"""Stores preallocated structs for use in Bloch CPU run_spin_precession function."""
mutable struct BlochCPUPrealloc{T} <: PreallocResult{T}
    B_old::VectorSU2{T}         # Vector{T}(Nspins x 1)
    B_new::VectorSU2{T}         # Vector{T}(Nspins x 1)
    φ::AbstractVector{T}        # Vector{T}(Nspins x 1)
    Rot::Spinor{T}              # Spinor{T}
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    BlochCPUPrealloc(
        @view(p.B_old[i]),
        @view(p.B_new[i]),
        @view(p.φ[i]),
        @view(p.Rot[i])
    )
end

"""Preallocates arrays for use in run_spin_precession."""
function prealloc(sim_method::Bloch, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}) where {T<:Real}
    return BlochCPUPrealloc(
        VectorSU2(zero(M.xy), zero(M.z)), #TODO: change this to not use M.xy
        VectorSU2(zero(M.xy), zero(M.z)),
        zero(obj.x),
        Spinor(zero(M.xy), zero(M.xy))
    )
end

"""
    run_spin_precession(obj, seq, Xt, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_precession! function in BlochSimpleSimulationMethod.jl 
optimized for the CPU. Uses a loop to step through time instead of allocating a matrix of size 
NSpins x seq.t. The Bz_old, Bz_new, ϕ, and Mxy arrays are pre-allocated in run_sim_time_iter! so 
that they can be re-used from block to block.
"""
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.CPU,
    pre::BlochCPUPrealloc
) where {T<:Real}
    # Initialize
    ADC_idx = 1
    Δt = zero(T) # Time
    pre.φ .= zero(T)
    # Simulation    
    for i in 1:length(seq.Δt)
        # Bz_eff(tn+1)
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[i+1,:]')
        @. pre.B_new.z = x * seq.Gx[i+1] + y * seq.Gy[i+1] + z * seq.Gz[i+1] + p.Δw / T(2π * γ)
        # Rotation angle
        @. pre.φ += T(-2π * γ) * (pre.B_old.z + pre.B_new.z) * seq.Δt[i] / 2
        Δt += seq.Δt[i]
        # Acquire signal
        if seq.ADC[i+1] # TODO: || sim_method.sample_all
            # Decay + Rot
            @. M.xy *= exp(-Δt / p.T2) * cis(pre.φ)
            # Sample
            sig[ADC_idx] = sum(M.xy) 
            ADC_idx += 1
            # Reset, decay and rot operators
            Δt = zero(T)
            pre.φ .= zero(T)
        end
        # Update
        pre.B_old.z .= pre.B_new.z
    end
    pre.B_old.xy .= zero(Complex{T})
    # Final Spin-State
    t_total = sum(seq.Δt)
    @. M.xy *= exp(-Δt / p.T2) * cis(pre.φ)
    @. M.z = M.z * exp(-t_total / p.T1) + p.ρ * (1 - exp(-t_total / p.T1))
    return nothing
end

"""
    run_spin_excitation!(obj, seq, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_excitation! function in BlochSimpleSimulationMethod.jl 
optimized for the CPU. Uses preallocation for all arrays to reduce memory usage.
"""
function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.CPU,
    pre::BlochCPUPrealloc
) where {T<:Real}
    # Init
    sample = 1
    # Simulation
    for i in 1:length(seq.Δt)
        # Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[i+1, :])
        # Effective field
        @. pre.B_new.xy = seq.B1[i+1]
        @. pre.B_new.z  = (seq.Gx[i+1] * x + seq.Gy[i+1] * y + seq.Gz[i+1] * z) + p.Δw / T(2π * γ) - seq.Δf[i+1] / T(γ)
        # Rotation
        calculateRot!(pre, seq.Δt[i])
        mul!(pre.Rot, M)
        # Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z  = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (1 - exp(-seq.Δt[i] / p.T1))
        # Sample
        if seq.ADC[i+1]
            sig[sample] .= sum(M.xy)
            sample += 1
        end
        # Update
        pre.B_old.xy .= pre.B_new.xy
        pre.B_old.z  .= pre.B_new.z
    end
    return nothing
end