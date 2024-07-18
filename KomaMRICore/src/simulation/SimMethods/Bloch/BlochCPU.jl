"""Stores arrays for use in Bloch CPU run_spin_precession function."""
struct BlochCPUPrealloc{T} <: PreallocResult{T}
    Bz_1::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bz_2::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bz_3::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bz_4::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bz_5::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Mxy_1::AbstractVector{Complex{T}}     # Vector{Complex{T}}(Nspins x 1)
    Mxy_2::AbstractVector{Complex{T}}     # Vector{Complex{T}}(Nspins x 1)
    Mxy_3::AbstractVector{Complex{T}}     # Vector{Complex{T}}(Nspins x 1)
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.Bz_1[i],
        p.Bz_2[i],
        p.Bz_3[i],
        p.Bz_4[i],
        p.Bz_5[i],
        p.Mxy_1[i],
        p.Mxy_2[i],
        p.Mxy_3[i]
    )
end

"""Preallocates arrays for use in run_spin_precession."""
function prealloc(sim_method::Bloch, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}) where {T<:Real}
    return BlochCPUPrealloc(
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
        similar(M.xy),
        similar(M.xy),
        similar(M.xy)
    )
end

"""
    run_spin_precession(obj, seq, Xt, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_precession! function in SimulationMethod.jl optimized
for the CPU. Uses a loop to step through time instead of allocating a matrix of size NSpins x seq.t.
The Bz_old, Bz_new, ϕ, and Mxy arrays are pre-allocated in run_sim_time_iter! so that they can be
re-used from block to block.
"""
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.CPU,
    prealloc::BlochCPUPrealloc
) where {T<:Real}
    #Simulation
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[1,:]')
    
    #Initialize arrays
    Bz_old = prealloc.Bz_1
    Bz_new = prealloc.Bz_2
    ϕ = prealloc.Bz_3
    Mxy = prealloc.Mxy_1
    fill!(ϕ, zero(T))
    @. Bz_old = x[:,1] * seq.Gx[1] + y[:,1] * seq.Gy[1] + z[:,1] * seq.Gz[1] + p.Δw / T(2π * γ)

    # Fill sig[1] if needed
    ADC_idx = 1
    if (seq.ADC[1])
        sig[1] = sum(M.xy)
        ADC_idx += 1
    end

    t_seq = zero(T) # Time
    for seq_idx=2:length(seq.t)
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[seq_idx,:]')
        t_seq += seq.Δt[seq_idx-1]

        #Effective Field
        @. Bz_new = x * seq.Gx[seq_idx] + y * seq.Gy[seq_idx] + z * seq.Gz[seq_idx] + p.Δw / T(2π * γ)
        
        #Rotation
        @. ϕ += (Bz_old + Bz_new) * T(-π * γ) * seq.Δt[seq_idx-1]

        #Acquired Signal
        if seq_idx <= length(seq.ADC) && seq.ADC[seq_idx]
            @. Mxy = exp(-t_seq / p.T2) * M.xy * cis(ϕ)
            sig[ADC_idx] = sum(Mxy) 
            ADC_idx += 1
        end

        Bz_old, Bz_new = Bz_new, Bz_old
    end

    #Final Spin-State
    @. M.xy = M.xy * exp(-t_seq / p.T2) * cis(ϕ)
    @. M.z = M.z * exp(-t_seq / p.T1) + p.ρ * (T(1) - exp(-t_seq / p.T1))

    return nothing
end

"""
    run_spin_excitation!(obj, seq, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_excitation! function in SimulationMethod.jl 
optimized for the CPU. Uses preallocation for all arrays to reduce memory usage.
"""
function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.CPU,
    prealloc::BlochCPUPrealloc
) where {T<:Real}
    ΔBz = prealloc.Bz_1
    Bz = prealloc.Bz_2
    B = prealloc.Bz_3
    φ = prealloc.Bz_4
    α = prealloc.Mxy_1
    β = prealloc.Mxy_2
    Maux_xy = prealloc.Mxy_3
    Maux_z = prealloc.Bz_5

    #Can be calculated outside of loop
    @. ΔBz = p.Δw / T(2π * γ)

    #Simulation
    for s in seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, s.t)
        #Effective field
        @. Bz = (s.Gx * x + s.Gy * y + s.Gz * z) + ΔBz - s.Δf / T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        @. B = sqrt(abs(s.B1)^2 + abs(Bz)^2)
        @. B[B == 0] = eps(T)
        #Spinor Rotation
        @. φ = T(-π * γ) * (B * s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler 
        @. α = cos(φ) - Complex{T}(im) * (Bz / B) * sin(φ)
        @. β = -Complex{T}(im) * (s.B1 / B) * sin(φ)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        @. M.xy = M.xy * exp(-s.Δt / p.T2)
        @. M.z = M.z * exp(-s.Δt / p.T1) + p.ρ * (T(1) - exp(-s.Δt / p.T1))
    end
    #Acquired signal
    #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end