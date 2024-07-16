"""Stores arrays for use in Bloch CPU run_spin_precession function."""
struct BlochCPUPrealloc{T} <: PreallocResult{T}
    Bz_old::AbstractVector{T}
    Bz_new::AbstractVector{T}
    ϕ::AbstractVector{T}
    Mxy::AbstractVector{Complex{T}}
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Mxy[i]
    )
end

"""Preallocates arrays for use in run_spin_precession."""
function prealloc(sim_method::Bloch, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}) where {T<:Real}
    return BlochCPUPrealloc(
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
        similar(M.xy)
    )
end

"""
    run_spin_precession(obj, seq, Xt, sig)

Specialized run_spin_precession function for simulating a sequence on the CPU with Bloch as the
simulation method. Compared with the default run_spin_precession function in SimMethod.jl, this
function uses a loop to step through time and does not allocate a matrix of size NSpins x seq.t.
The four arrays of size Nspins x 1 are preallocated in SimulatorCore::run_sim_time_iter! so that
they can be re-used from block-to-block.
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
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    
    #Initialize arrays
    Bz_old = prealloc.Bz_old
    Bz_new = prealloc.Bz_new
    ϕ = prealloc.ϕ
    Mxy = prealloc.Mxy
    fill!(ϕ, zero(T))
    Bz_old .= x[:,1] .* seq.Gx[1] .+ y[:,1] .* seq.Gy[1] .+ z[:,1] .* seq.Gz[1] .+ p.Δw / T(2π * γ)

    # Fill sig[1] if needed
    ADC_idx = 1
    if (seq.ADC[1])
        sig[1] = sum(M.xy)
        ADC_idx += 1
    end

    t_seq = zero(T) # Time
    for seq_idx=2:length(seq.t)
        t_seq += seq.Δt[seq_idx-1]

        #Effective Field
        if size(x,2) > 1 #Motion
            Bz_new .= x[:,seq_idx] .* seq.Gx[seq_idx] .+ y[:,seq_idx] .* seq.Gy[seq_idx] .+ z[:,seq_idx] .* seq.Gz[seq_idx] .+ p.Δw / T(2π * γ)
        else             #No motion
            Bz_new .= x .* seq.Gx[seq_idx] .+ y .* seq.Gy[seq_idx] .+ z.* seq.Gz[seq_idx] .+ p.Δw / T(2π * γ)
        end
        
        #Rotation
        ϕ .= ϕ .+ (Bz_old .+ Bz_new) .* (T(-2π * γ) * seq.Δt[seq_idx-1] / 2)

        #Acquired Signal
        if seq_idx <= length(seq.ADC) && seq.ADC[seq_idx]
            Mxy .= exp.(-t_seq ./ p.T2) .* (M.xy .* (cos.(ϕ) .+ im * sin.(ϕ)))
            sig[ADC_idx] = sum(Mxy) 
            ADC_idx += 1
        end

        Bz_old, Bz_new = Bz_new, Bz_old
    end

    #Final Spin-State
    M.xy .= M.xy .* exp.(-t_seq ./ p.T2) .* (cos.(ϕ) .+ im * sin.(ϕ))
    M.z .= M.z .* exp.(-t_seq ./ p.T1) .+ p.ρ .* (1 .- exp.(-t_seq ./ p.T1))

    return nothing
end
