include("KernelFunctions.jl")
include("PrecessionKernel.jl")
include("ExcitationKernel.jl")

"""Stores preallocated arrays for use in Bloch GPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochGPUPrealloc{T} <: PreallocResult{T}
    sig_output::AbstractMatrix{Complex{T}}
    ΔBz::AbstractVector{T}
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(
    sim_method::SM, 
    backend::KA.GPU, 
    obj::Phantom{T}, 
    M::Mag{T}, 
    max_block_length::Integer, 
    groupsize
) where {T<:Real, SM<:BlochLikeSimMethods}
    return BlochGPUPrealloc(
        KA.zeros(backend, Complex{T}, (cld(size(obj.x, 1), groupsize), max_block_length)),
        obj.Δw ./ T(2π .* γ)
    )
end

prealloc(
    sim_method::BlochMagnusBGL4,
    backend::KA.GPU,
    obj::Phantom{T},
    M::Mag{T},
    max_block_length::Integer,
    groupsize
) where {T<:Real} =
    BlochGPUPrealloc(
        KA.zeros(backend, Complex{T}, (cld(size(obj.x, 1), groupsize), max_block_length)),
        obj.Δw ./ T(2π .* γ)
    )

prealloc(
    sim_method::BlochMagnusBGL6,
    backend::KA.GPU,
    obj::Phantom{T},
    M::Mag{T},
    max_block_length::Integer,
    groupsize
) where {T<:Real} =
    prealloc(BlochMagnusBGL4(), backend, obj, M, max_block_length, groupsize)

function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
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
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), reduction_mode(backend), Val(has_adc),
        BlochMagnusConst1(),
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        reduce_signal_output!(sig, view(pre.sig_output, :, 1:length(sig)), backend)
    end

    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
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
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), reduction_mode(backend), Val(has_adc),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    if has_adc
        reduce_signal_output!(sig, view(pre.sig_output, :, 1:length(sig)), backend)
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
    sim_method::BlochMagnusBGL6,
    groupsize::Integer,
    backend::KA.Backend,
    pre::BlochGPUPrealloc
) where {T<:Real} =
    run_spin_precession!(p, seq, sig, M, BlochMagnusBGL4(), groupsize, backend, pre)

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
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
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), reduction_mode(backend), Val(has_adc),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        reduce_signal_output!(sig, view(pre.sig_output, :, 1:length(sig)), backend)
    end

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
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
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), reduction_mode(backend), Val(has_adc),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    if has_adc
        reduce_signal_output!(sig, view(pre.sig_output, :, 1:length(sig)), backend)
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
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ψ, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), reduction_mode(backend), Val(has_adc),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    if has_adc
        reduce_signal_output!(sig, view(pre.sig_output, :, 1:length(sig)), backend)
    end

    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end
