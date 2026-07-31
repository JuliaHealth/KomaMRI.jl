using Test
using Enzyme
using FiniteDifferences
using KomaMRIBase
using KomaMRICore
using Reactant

Reactant.set_default_backend("cuda")
Reactant.allowscalar(false)

const REACTANT_CUDA_RF0 = [1.3, 1.7, 1.1]

function reactant_cuda_traced_vector(template, values, ::Type{T}=Float64) where {T}
    out = similar(template, T, length(values))
    Reactant.allowscalar() do
        for i in eachindex(values)
            out[i] = values[i]
        end
    end
    return out
end

function reactant_cuda_discrete_sequence(rf, t_start, t_end; excitation)
    n = length(rf)
    t = collect(range(t_start, t_end; length=n))
    z = zero.(rf)
    B1 = excitation ? complex.(rf) .* 1e-6 : complex.(z)
    return DiscreteSequence(
        z,
        copy(z),
        copy(z),
        B1,
        copy(z),
        copy(z),
        falses(n),
        fill(excitation, n - 1),
        t,
        diff(t),
    )
end

function reactant_cuda_loss(rf)
    obj = Phantom(
        x=reactant_cuda_traced_vector(rf, [0.0, 1e-2]),
        ρ=reactant_cuda_traced_vector(rf, [1.0, 0.75]),
        T1=reactant_cuda_traced_vector(rf, [1.0, 0.8]),
        T2=reactant_cuda_traced_vector(rf, [0.08, 0.12]),
        Δw=reactant_cuda_traced_vector(rf, 2π .* [10.0, -12.0]),
    )
    M = Mag(
        reactant_cuda_traced_vector(rf, zeros(2), ComplexF64),
        reactant_cuda_traced_vector(rf, [1.0, 0.75]),
    )
    sig = similar(rf, ComplexF64, 0, 1, 1)
    method = BlochSimple()
    # KernelAbstractions selects the traceable Koma path; Reactant/XLA selects CUDA.
    backend = KomaMRICore.KA.CPU()
    prealloc = KomaMRICore.DefaultPrealloc()

    KomaMRICore.run_spin_excitation_parallel!(
        obj,
        reactant_cuda_discrete_sequence(rf, 0.0, 0.6e-3; excitation=true),
        sig,
        M,
        method,
        256,
        backend,
        prealloc;
        Nthreads=1,
    )
    KomaMRICore.run_spin_precession_parallel!(
        obj,
        reactant_cuda_discrete_sequence(rf, 0.6e-3, 1.4e-3; excitation=false),
        sig,
        M,
        method,
        256,
        backend,
        prealloc;
        Nthreads=1,
    )

    target = reactant_cuda_traced_vector(rf, [0.4, 0.2])
    return sum(abs2, M.xy) + sum(abs2, M.z .- target)
end

reactant_cuda_gradient(rf) =
    Enzyme.gradient(Enzyme.ReverseWithPrimal, reactant_cuda_loss, rf).derivs[1]

function reactant_cuda_sequence(rf_scale)
    rf_duration = 0.6e-3
    total_duration = 1.4e-3
    adc_delay = 0.0
    adc_duration = total_duration - adc_delay
    rf = RF(
        complex.(rf_scale) .* 1e-6,
        rf_duration,
        0.0,
        0.0,
        rf_duration / 2,
        0.0,
        Excitation(),
        Val(:preserve),
    )
    gradient = Grad(0.0, 0.0)
    adc = ADC(3, adc_duration, adc_delay)
    return Sequence(
        reshape([gradient, gradient, gradient], 3, 1),
        reshape([rf], 1, 1),
        [adc],
        [total_duration],
        [Extension[]],
        Dict{String,Any}(),
    )
end

function reactant_cuda_discretize_loss(rf_scale)
    seqd = discretize(
        reactant_cuda_sequence(rf_scale);
        sampling_rule=MaxStepSizeRule(1e-3, 0.2e-3),
    )
    parts, excitation_bool = KomaMRICore.get_sim_ranges(seqd)
    zeros2 = zero.(rf_scale[1:2])
    density = reactant_cuda_traced_vector(rf_scale, [1.0, 0.75])
    obj = Phantom(
        x=copy(zeros2),
        ρ=density,
        T1=one.(density),
        T2=0.1 .* one.(density),
        Δw=copy(zeros2),
    )
    M = Mag(complex.(zeros2), copy(density))
    signal = similar(rf_scale, ComplexF64, sum(seqd.ADC), 1, 1)
    KomaMRICore.run_sim_time_iter!(
        obj,
        seqd,
        signal,
        M,
        BlochSimple(),
        KomaMRICore.KA.CPU();
        parts,
        excitation_bool,
        Nblocks=length(parts),
        Nthreads=1,
    )
    return sum(abs2, signal)
end

reactant_cuda_discretize_gradient(rf) =
    Enzyme.gradient(Enzyme.ReverseWithPrimal, reactant_cuda_discretize_loss, rf).derivs[1]

@testset "Reactant + Enzyme CUDA through parallel BlochSimple kernels" begin
    devices = Reactant.devices()
    @info "Reactant CUDA devices" devices=string.(devices)
    @test !isempty(devices)
    @test all(occursin("CUDA", uppercase(string(device))) for device in devices)

    rf = Reactant.to_rarray(REACTANT_CUDA_RF0)
    compiled_loss = Reactant.@compile sync=true reactant_cuda_loss(rf)
    compiled_gradient = Reactant.@compile sync=true reactant_cuda_gradient(rf)
    finite_difference = FiniteDifferences.grad(
        FiniteDifferences.central_fdm(5, 1), reactant_cuda_loss, REACTANT_CUDA_RF0,
    )[1]

    @test Reactant.to_number(compiled_loss(rf)) ≈ reactant_cuda_loss(REACTANT_CUDA_RF0)
    @test Array(compiled_gradient(rf)) ≈ finite_difference rtol=1e-8 atol=1e-10
end

@testset "Reactant + Enzyme CUDA through discretize, ADC, and BlochSimple iterator" begin
    rf = Reactant.to_rarray(REACTANT_CUDA_RF0)
    compiled_loss = Reactant.@compile sync=true reactant_cuda_discretize_loss(rf)
    compiled_gradient = Reactant.@compile sync=true reactant_cuda_discretize_gradient(rf)
    finite_difference = FiniteDifferences.grad(
        FiniteDifferences.central_fdm(5, 1),
        reactant_cuda_discretize_loss,
        REACTANT_CUDA_RF0,
    )[1]

    @test Reactant.to_number(compiled_loss(rf)) ≈
        reactant_cuda_discretize_loss(REACTANT_CUDA_RF0)
    @test Array(compiled_gradient(rf)) ≈ finite_difference rtol=1e-8 atol=1e-10
end
