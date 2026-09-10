include("KernelFunctions.jl")
include("PrecessionKernel.jl")
include("ExcitationKernel.jl")

"""Stores preallocated arrays for use in Bloch GPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochGPUPrealloc{
    SignalType<:AbstractMatrix,
    ReducedSignalType<:AbstractMatrix,
    OffResonanceType<:AbstractVector,
    S,
    P,
} <: PreallocResult
    sig_output::SignalType
    sig_output_final::ReducedSignalType
    ΔBz::OffResonanceType
    receiver::S
    coordinates::P
end

function bloch_gpu_prealloc(
    backend, obj, max_block_length, max_adc_samples, groupsize, sys,
)
    T = eltype(obj.ρ)
    ncoils = get_n_coils(sys.receiver)
    signal_length = max_adc_samples * ncoils
    return BlochGPUPrealloc(
        KA.zeros(backend, Complex{T}, cld(length(obj), groupsize), signal_length),
        KA.zeros(backend, Complex{T}, 1, signal_length),
        obj.Δw ./ T(2π .* γ),
        prealloc_receiver(sys.receiver, obj, backend, obj.motion),
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(
    ::BlochLikeSimMethods,
    backend::KA.GPU,
    obj,
    _M,
    max_block_length,
    max_adc_samples,
    groupsize,
    sys,
)
    return bloch_gpu_prealloc(
        backend, obj, max_block_length, max_adc_samples, groupsize, sys,
    )
end

prealloc(
    ::BlochMagnusBGL4,
    backend::KA.GPU,
    obj,
    _M,
    max_block_length,
    max_adc_samples,
    groupsize,
    sys,
) =
    bloch_gpu_prealloc(
        backend, obj, max_block_length, max_adc_samples, groupsize, sys,
    )

prealloc(
    ::BlochMagnusBGL6,
    backend::KA.GPU,
    obj,
    M,
    max_block_length,
    max_adc_samples,
    groupsize,
    sys,
) =
    prealloc(
        BlochMagnusBGL4(), backend, obj, M, max_block_length, max_adc_samples,
        groupsize, sys,
    )

function reduce_signal_groups!(sig, pre, nspins, groupsize)
    groups = 1:cld(nspins, groupsize)
    samples = 1:length(sig)
    AK.reduce(
        +, view(pre.sig_output, groups, samples);
        init=zero(eltype(sig)), dims=1,
        temp=view(pre.sig_output_final, :, samples),
    )
    sig .= reshape(view(pre.sig_output_final, 1, samples), size(sig))
    return nothing
end

function run_spin_precession!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::BlochMagnusBGL4,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
)
    T = eltype(p.ρ)
    x, y, z = spin_coordinates!(
        pre.coordinates, p.motion, p.x, p.y, p.z, seq.t',
    )
    has_adc = !isempty(sig)

    precession_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.receiver, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        motion_enabled(pre.coordinates), Val(supports_warp_reduction(backend)),
        Val(has_adc), has_coil_sensitivities(pre.receiver),
        BlochMagnusConst1(),
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    has_adc && reduce_signal_groups!(sig, pre, length(M.xy), groupsize)

    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_precession!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::SM,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {SM<:BlochLikeSimMethods}
    T = eltype(p.ρ)
    #Motion
    x, y, z = spin_coordinates!(
        pre.coordinates, p.motion, p.x, p.y, p.z, seq.t',
    )
    has_adc = !isempty(sig)

    #Precession
    precession_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.receiver, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        motion_enabled(pre.coordinates), Val(supports_warp_reduction(backend)),
        Val(has_adc), has_coil_sensitivities(pre.receiver),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    has_adc && reduce_signal_groups!(sig, pre, length(M.xy), groupsize)

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)

    return nothing
end

run_spin_precession!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::BlochMagnusBGL6,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) =
    run_spin_precession!(p, seq, sig, M, sys, BlochMagnusBGL4(), groupsize, backend, pre)

function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::BlochMagnusBGL4,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
)
    T = eltype(p.ρ)
    x, y, z = spin_coordinates!(
        pre.coordinates, p.motion, p.x, p.y, p.z, seq.t',
    )
    has_adc = !isempty(sig)

    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.receiver, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        motion_enabled(pre.coordinates), Val(supports_warp_reduction(backend)),
        Val(has_adc), has_coil_sensitivities(pre.receiver),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    has_adc && reduce_signal_groups!(sig, pre, length(M.xy), groupsize)

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::SM,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {SM<:BlochLikeSimMethods}
    T = eltype(p.ρ)
    #Motion
    x, y, z = spin_coordinates!(
        pre.coordinates, p.motion, p.x, p.y, p.z, seq.t',
    )
    has_adc = !isempty(sig)

    #Excitation
    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.receiver, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        motion_enabled(pre.coordinates), Val(supports_warp_reduction(backend)),
        Val(has_adc), has_coil_sensitivities(pre.receiver),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    has_adc && reduce_signal_groups!(sig, pre, length(M.xy), groupsize)

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ) # TODO: reset state inside kernel

    return nothing
end

run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::BlochMagnusBGL6,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) =
begin
    T = eltype(p.ρ)
    x, y, z = spin_coordinates!(
        pre.coordinates, p.motion, p.x, p.y, p.z, seq.t',
    )
    has_adc = !isempty(sig)

    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.receiver, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        motion_enabled(pre.coordinates), Val(supports_warp_reduction(backend)),
        Val(has_adc), has_coil_sensitivities(pre.receiver),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    has_adc && reduce_signal_groups!(sig, pre, length(M.xy), groupsize)

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end
