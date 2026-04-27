include("KernelFunctions.jl")
include("PrecessionKernel.jl")
include("ExcitationKernel.jl")

"""Stores preallocated arrays for use in Bloch GPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochGPUPrealloc{T} <: PreallocResult{T}
    sig_output::AbstractMatrix{Complex{T}}
    sig_output_final::AbstractMatrix{Complex{T}}
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
        KA.zeros(backend, Complex{T}, 1, max_block_length),
        obj.Δw ./ T(2π .* γ)
    )
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

    #Precession
    precession_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
    sig .= transpose(view(pre.sig_output_final,:,1:length(sig)))

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)

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

    #Excitation
    excitation_kernel!(backend, groupsize)(
        pre.sig_output,
        M.xy, M.z,
        x, y, z, pre.ΔBz, p.T1, p.T2, p.ρ, UInt32(length(M.xy)),
        seq.Gx, seq.Gy, seq.Gz, seq.Δt, seq.Δf, seq.B1, seq.ADC, UInt32(length(seq.t)),
        Val(!(p.motion isa NoMotion)), Val(supports_warp_reduction(backend)),
        sim_method,
        ndrange=(cld(length(M.xy), groupsize) * groupsize)
    )

    #Signal
    AK.reduce(+, view(pre.sig_output,:,1:length(sig)); init=zero(Complex{T}), dims=1, temp=view(pre.sig_output_final,:,1:length(sig)))
    sig .= transpose(view(pre.sig_output_final,:,1:length(sig)))

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ) # TODO: reset state inside kernel

    return nothing
end