using FiniteDifferences: central_fdm, grad

const BLOCHSIMPLE_AD_RF0 = [1.3, 1.7, 1.1]
const BLOCHSIMPLE_AD_DIRECTION = [0.2, -0.1, 0.15]

function blochsimple_ad_sequence(rf_scale)
    Trf = 0.6e-3
    Tadc = 0.8e-3
    seq = Sequence()
    @addblock seq += RF(complex.(rf_scale) .* 1e-6, Trf)
    @addblock seq += ADC(4, Tadc)
    return seq
end

function blochsimple_ad_loss(rf_scale)
    obj = Phantom(
        x=[0.0, 1e-2],
        ρ=[1.0, 0.75],
        T1=[1.0, 0.8],
        T2=[0.08, 0.12],
        Δw=2π .* [10.0, -12.0],
    )
    sim_params = Dict{String, Any}(
        "sim_method" => KomaMRICore.BlochSimple(),
        "gpu" => false,
        "Nthreads" => 1,
        "return_type" => "state",
        "precision" => "f64",
        "Δt_rf" => 0.2e-3,
    )
    M = simulate(obj, blochsimple_ad_sequence(rf_scale), Scanner(); sim_params, verbose=false)
    target_z = [0.4, 0.2]
    return sum(abs2, M.xy) + sum(abs2, M.z .- target_z)
end

blochsimple_ad_fd_gradient(rf_scale=BLOCHSIMPLE_AD_RF0) =
    grad(central_fdm(5, 1), blochsimple_ad_loss, rf_scale)[1]

function blochsimple_ad_gradient_matches_fd(ad_grad; rf_scale=BLOCHSIMPLE_AD_RF0)
    fd_grad = blochsimple_ad_fd_gradient(rf_scale)
    return all(isfinite, ad_grad) && isapprox(ad_grad, fd_grad; rtol=1e-3, atol=1e-7)
end

const BLOCHSIMPLE_DISCRETIZE_AD_RF0 = [1.3, 1.7, 1.1]

function blochsimple_discretize_ad_vector(template, values, ::Type{T}=Float64) where {T}
    out = similar(template, T, length(values))
    Reactant.allowscalar() do
        for i in eachindex(values)
            out[i] = values[i]
        end
    end
    return out
end

function blochsimple_discretize_ad_sequence(rf_scale; adc_delay=0.0)
    rf_duration = 0.6e-3
    total_duration = 1.4e-3
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

function blochsimple_discretize_ad_loss(rf_scale)
    seqd = discretize(
        blochsimple_discretize_ad_sequence(rf_scale);
        sampling_rule=MaxStepSizeRule(1e-3, 0.2e-3),
    )
    parts, excitation_bool = KomaMRICore.get_sim_ranges(seqd)
    zeros2 = zero.(rf_scale[1:2])
    density = blochsimple_discretize_ad_vector(rf_scale, (1.0, 0.75))
    obj = Phantom(
        x=copy(zeros2),
        ρ=density,
        T1=one.(density),
        T2=0.1 .* one.(density),
        Δw=copy(zeros2),
    )
    magnetization = KomaMRICore.Mag(complex.(zeros2), copy(density))
    signal = similar(rf_scale, ComplexF64, sum(seqd.ADC), 1, 1)
    KomaMRICore.run_sim_time_iter!(
        obj,
        seqd,
        signal,
        magnetization,
        KomaMRICore.BlochSimple(),
        KomaMRICore.KA.CPU();
        parts,
        excitation_bool,
        Nblocks=length(parts),
        Nthreads=1,
    )
    return sum(abs2, signal)
end

blochsimple_discretize_ad_fd_gradient(rf_scale=BLOCHSIMPLE_DISCRETIZE_AD_RF0) =
    grad(central_fdm(5, 1), blochsimple_discretize_ad_loss, rf_scale)[1]
