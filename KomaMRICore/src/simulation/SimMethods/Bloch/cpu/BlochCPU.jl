"""Stores preallocated structs for use in Bloch CPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochCPUPrealloc{T} <: PreallocResult{T}
    M::Mag{T}                               # Mag{T}
    Bz_old::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bz_new::AbstractVector{T}               # Vector{T}(Nspins x 1)
    ϕ::AbstractVector{T}                    # Vector{T}(Nspins x 1)
    Rot::Spinor{T}                          # Spinor{T}
    ΔBz::AbstractVector{T}            # Vector{T}(Nspins x 1)
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.M[i],
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Rot[i],
        p.ΔBz[i]
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
        obj.Δw ./ T(2π .* γ)
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
    seq::AbstractDiscreteSequence,
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
    #Initialize
    fill!(ϕ, zero(T))
    block_time = zero(T)
    sample = 1
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        #Effective Field
        get_Bz_field!(Bz_new, seq, x, y, z, i + 1)
        Bz_new .+= ΔBz
        #Rotation
        @. ϕ += (Bz_old + Bz_new) * T(-π * γ) * seq.Δt[i]
        block_time += seq.Δt[i]
        #Acquired Signal
        if seq.ADC[i + 1]
            #Update signal
            @. Mxy = exp(-block_time / p.T2) * M.xy * cis(ϕ)
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
    @. M.xy = M.xy * exp(-block_time / p.T2) * cis(ϕ)
    @. M.z = M.z * exp(-block_time / p.T1) + p.ρ * (T(1) - exp(-block_time / p.T1))
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
    seq::AbstractDiscreteSequence{T},
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
    #Initialize
    sample = 1
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[i])
        #Effective field
        get_Bz_field!(Bz, seq, x, y, z, i)
        @. Bz = Bz + ΔBz - seq.Δf[i] / T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        @. B = sqrt(abs(seq.B1[i])^2 + abs(Bz)^2)
        @. B[B == 0] = eps(T)
        #Spinor Rotation
        @. φ_half = T(-π * γ) * (B * seq.Δt[i]) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        @. α = cos(φ_half) - Complex{T}(im) * (Bz / B) * sin(φ_half)
        @. β = -Complex{T}(im) * (seq.B1[i] / B) * sin(φ_half)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset!(M,  seq.t[i + 1, :], p.motion; replace_by=p.ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            sig[sample] = sum(M.xy) 
            sample += 1
        end
    end
    return nothing
end