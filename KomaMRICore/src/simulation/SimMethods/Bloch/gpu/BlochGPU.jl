include("KernelFunctions.jl")
include("PrecessionKernel.jl")
include("ExcitationKernel.jl")

"""Stores preallocated arrays for use in Bloch GPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochGPUPrealloc{T,S,H} <: PreallocResult{T}
    sig_output::AbstractMatrix{Complex{T}}
    sig_output_final::AbstractMatrix{Complex{T}}
    ΔBz::AbstractVector{T}
    sens::S
    has_sens::H
end

gpu_sens(::UniformCoilSens, obj, backend) = nothing, Val(false)
gpu_sens(receiver::BirdcageCoilSens, obj, backend) =
    get_sens(receiver, obj.x, obj.y, obj.z), Val(true)
gpu_sens(receiver::ArbitraryCoilSens, obj, backend) =
    gpu(get_sens(cpu(receiver), cpu(obj.x), cpu(obj.y), cpu(obj.z)), backend), Val(true)

function bloch_gpu_prealloc(backend, obj::Phantom{T}, max_block_length, groupsize, sys) where {T}
    ncoils = get_n_coils(sys.receiver)
    sens, has_sens = gpu_sens(sys.receiver, obj, backend)
    signal_length = max_block_length * ncoils
    return BlochGPUPrealloc(
        KA.zeros(backend, Complex{T}, cld(length(obj), groupsize), signal_length),
        KA.zeros(backend, Complex{T}, 1, signal_length),
        obj.Δw ./ T(2π .* γ),
        sens,
        has_sens,
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(
    sim_method::SM, 
    backend::KA.GPU, 
    obj::Phantom{T}, 
    M::Mag{T}, 
    max_block_length::Integer, 
    groupsize,
    sys::Scanner,
) where {T<:Real, SM<:BlochLikeSimMethods}
    return bloch_gpu_prealloc(backend, obj, max_block_length, groupsize, sys)
end

prealloc(
    sim_method::BlochMagnusBGL4,
    backend::KA.GPU,
    obj::Phantom{T},
    M::Mag{T},
    max_block_length::Integer,
    groupsize,
    sys::Scanner,
) where {T<:Real} =
    bloch_gpu_prealloc(backend, obj, max_block_length, groupsize, sys)

prealloc(
    sim_method::BlochMagnusBGL6,
    backend::KA.GPU,
    obj::Phantom{T},
    M::Mag{T},
    max_block_length::Integer,
    groupsize,
    sys::Scanner,
) where {T<:Real} =
    prealloc(BlochMagnusBGL4(), backend, obj, M, max_block_length, groupsize, sys)

function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnusBGL4,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real}
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    has_adc = length(sig) > 0

    precession_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.sens, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)), Val(has_adc), pre.has_sens,
        BlochMagnusConst1(),
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
        sig .= reshape(view(pre.sig_output_final, 1, 1:length(sig)), size(sig))
    end

    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::SM,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real, SM<:BlochLikeSimMethods}
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    has_adc = length(sig) > 0

    #Precession
    precession_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.sens, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)), Val(has_adc), pre.has_sens,
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    if has_adc
        AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
        sig .= reshape(view(pre.sig_output_final, 1, 1:length(sig)), size(sig))
    end

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)

    return nothing
end

run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnusBGL6,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real} =
    run_spin_precession!(p, seq, sig, M, sys, BlochMagnusBGL4(), groupsize, backend, pre)

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnusBGL4,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real}
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    has_adc = length(sig) > 0

    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.sens, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)), Val(has_adc), pre.has_sens,
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
        sig .= reshape(view(pre.sig_output_final, 1, 1:length(sig)), size(sig))
    end

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::SM,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real, SM<:BlochLikeSimMethods}
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    has_adc = length(sig) > 0

    #Excitation
    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.sens, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)), Val(has_adc), pre.has_sens,
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    if has_adc
        AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
        sig .= reshape(view(pre.sig_output_final, 1, 1:length(sig)), size(sig))
    end

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ) # TODO: reset state inside kernel

    return nothing
end

run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnusBGL6,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real} =
begin
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    has_adc = length(sig) > 0

    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        pre.sens, UInt32(size(sig, 2)), UInt32(size(sig, 1)),
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)), Val(has_adc), pre.has_sens,
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
        sig .= reshape(view(pre.sig_output_final, 1, 1:length(sig)), size(sig))
    end

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end
