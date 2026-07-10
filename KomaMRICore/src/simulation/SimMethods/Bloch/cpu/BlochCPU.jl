struct RelaxationCPUPrealloc{T,RV<:AbstractVector{T}}
    neg_inv_T1::RV
    neg_inv_T2::RV
    E1::RV
end

Base.view(p::RelaxationCPUPrealloc, i::UnitRange) =
    RelaxationCPUPrealloc(view(p.neg_inv_T1, i), view(p.neg_inv_T2, i), view(p.E1, i))

relaxation_prealloc(obj::Phantom) =
    RelaxationCPUPrealloc(.-inv.(obj.T1), .-inv.(obj.T2), similar(obj.T1))

"""Stores preallocated structs for use in Bloch CPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochCPUPrealloc{
    T,
    MT<:Mag{T},
    RV<:AbstractVector{T},
    ST<:Spinor{T},
} <: PreallocResult{T}
    M::MT
    Bz_old::RV
    Bz_new::RV
    ϕ::RV
    Rot::ST
    ΔBz::RV
    relaxation::RelaxationCPUPrealloc{T,RV}
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.M[i],
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Rot[i],
        p.ΔBz[i],
        view(p.relaxation, i)
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(sim_method::Bloch, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real}
    return BlochCPUPrealloc(
        Mag(
            similar(M.xy),
            similar(M.z)
        ),
        zeros(T, size(obj.x)),
        zeros(T, size(obj.x)),
        zeros(T, size(obj.x)),
        Spinor(
            similar(M.xy),
            similar(M.xy)
        ),
        obj.Δw ./ T(2π .* γ),
        relaxation_prealloc(obj)
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
    groupsize,
    backend::KA.CPU,
    prealloc::PreallocResult{T}
) where {T<:Real}
    #Rename arrays
    Bz_old = prealloc.Bz_old
    Bz_new = prealloc.Bz_new
    ϕ = prealloc.ϕ
    Mxy = prealloc.M.xy
    ΔBz = prealloc.ΔBz
    (; neg_inv_T1, neg_inv_T2, E1) = prealloc.relaxation
    #Initialize
    B_to_ω_half = T(-π * γ)
    fill!(ϕ, zero(T))
    block_time = zero(T)
    sample = 1
    x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[1])
    @. Bz_old = x * seq.Gx[1] + y * seq.Gy[1] + z * seq.Gz[1] + ΔBz
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        #Effective Field
        @. Bz_new = x * seq.Gx[i + 1] + y * seq.Gy[i + 1] + z * seq.Gz[i + 1] + ΔBz
        #Rotation
        @. ϕ += (Bz_old + Bz_new) * B_to_ω_half * seq.Δt[i]
        block_time += seq.Δt[i]
        #Acquired Signal
        if seq.ADC[i + 1]
            #Update signal
            @. Mxy = exp(block_time * neg_inv_T2) * M.xy * cis(ϕ)
            #Reset Spin-State (Magnetization). Only for FlowPath
            outflow_spin_reset!(Mxy, seq.t[i + 1], p.motion)
            #Acquired signal
            sig[sample] = sum(Mxy) 
            sample += 1
        end
        #Update simulation state
        Bz_old .= Bz_new
    end
    #Final Spin-State
    @. M.xy = M.xy * exp(block_time * neg_inv_T2) * cis(ϕ)
    @. E1 = exp(block_time * neg_inv_T1)
    @. M.z = M.z * E1 + p.ρ * (T(1) - E1)
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
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
    groupsize,
    backend::KA.CPU,
    prealloc::BlochCPUPrealloc
) where {T<:Real}
    #Rename arrays 
    Bz = prealloc.Bz_old
    B = prealloc.Bz_new
    φ_half = prealloc.ϕ
    α = prealloc.Rot.α
    β = prealloc.Rot.β
    ΔBz = prealloc.ΔBz
    Maux_xy = prealloc.M.xy
    Maux_z = prealloc.M.z
    (; neg_inv_T1, neg_inv_T2, E1) = prealloc.relaxation
    #Initialize
    B_to_ω_half = T(-π * γ)
    inv_γ = inv(T(γ))
    sample = 1
    # Rotating frame -> RF frame
    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i])
        #Effective field
        @. Bz = (seq.Gx[i] * x + seq.Gy[i] * y + seq.Gz[i] * z) + ΔBz - seq.Δf[i] * inv_γ # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        @. B = sqrt(abs2(seq.B1[i]) + Bz^2)
        #Spinor Rotation
        @. φ_half = B_to_ω_half * (B * seq.Δt[i]) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        @. α = cos(φ_half)
        @. B = sin(φ_half) / (B + (B == 0) * eps(T))
        @. α -= Complex{T}(im) * Bz * B
        @. β = -Complex{T}(im) * seq.B1[i] * B
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        @. M.xy = M.xy * exp(seq.Δt[i] * neg_inv_T2)
        @. E1 = exp(seq.Δt[i] * neg_inv_T1)
        @. M.z = M.z * E1 + p.ρ * (T(1) - E1)
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset_at!(M, seq.t, i + 1, p.motion; replace_by=p.ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            sig[sample] = sum(M.xy) 
            sample += 1
        end
    end
    # RF frame -> Rotating frame
    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
