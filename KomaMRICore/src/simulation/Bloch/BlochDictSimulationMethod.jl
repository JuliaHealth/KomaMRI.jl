Base.@kwdef struct BlochDict <: SimulationMethod
    save_Mz::Bool = false
end

export BlochDict
Base.show(io::IO, s::BlochDict) = begin
    print(io, "BlochDict(save_Mz=$(s.save_Mz))")
end

function sim_output_dim(
    obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochDict
) where {T<:Real}
    out_state_dim = sim_method.save_Mz ? 2 : 1
    return (sum(seq.ADC.N), length(obj), out_state_dim)
end

"""Preallocated arrays for use in run_spin_precession."""
struct BlochDictPrealloc{T} <: PreallocResult{T}
    Bz_1::AbstractVector{T}
    Bz_2::AbstractVector{T}
    ϕ::AbstractVector{T}
end

Base.view(p::BlochDictPrealloc, i::UnitRange) = begin
    @views BlochDictPrealloc(
        p.Bz_1[i],
        p.Bz_2[i],
        p.ϕ[i]
    )
end

"""BlochDict preallocation function. Returns arrays for use in run_spin_precession."""
function prealloc(sim_method::BlochDict, obj::Phantom{T}, M::Mag{T}) where {T<:Real}
    BlochDictPrealloc(
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
    )
end

"""
    run_spin_precession(obj, seq, Xt, sig)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
precession.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochDict,
    backend::KA.Backend,
    prealloc::BlochDictPrealloc
) where {T<:Real}
    #Simulation
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    
    #Initialize arrays
    Bz_1 = prealloc.Bz_1
    Bz_2 = prealloc.Bz_2
    ϕ = prealloc.ϕ
    fill!(ϕ, zero(T))
    Bz_1 .= x[:,1] .* seq.Gx[1] .+ y[:,1] .* seq.Gy[1] .+ z[:,1] .* seq.Gz[1] .+ p.Δw / T(2π * γ)

    # Fill sig[1] if needed
    ADC_idx = 1
    if (seq.ADC[1])
        sig[1,:,1] .= M.xy
        if sim_method.save_Mz
            sig[1,:,2] .= M.z
        end
        ADC_idx += 1
    end

    t_seq = zero(T) # Time
    for seq_idx=2:length(seq.t)
        t_seq += seq.Δt[seq_idx-1]

        #Effective Field
        if size(x,2) > 1 #Motion
            Bz_2 .= x[:,seq_idx] .* seq.Gx[seq_idx] .+ y[:,seq_idx] .* seq.Gy[seq_idx] .+ z[:,seq_idx] .* seq.Gz[seq_idx] .+ p.Δw / T(2π * γ)
        else             #No motion
            Bz_2 .= x .* seq.Gx[seq_idx] .+ y .* seq.Gy[seq_idx] .+ z.* seq.Gz[seq_idx] .+ p.Δw / T(2π * γ)
        end
        
        #Rotation
        ϕ .= ϕ .+ (Bz_1 .+ Bz_2) .* (T(-2π * γ) * seq.Δt[seq_idx-1] / 2)

        #Acquired Signal
        if seq_idx <= length(seq.ADC) && any(seq.ADC[seq_idx,:])
            sig[ADC_idx,:,1] .= M.xy .* exp.(-t_seq ./ p.T2) .* (cos.(ϕ) .+ im * sin.(ϕ))
            if sim_method.save_Mz
                sig[ADC_idx,:,2] .= M.z .* exp.(-t_seq ./ p.T1) .+ p.ρ .* (1 .- exp.(-t_seq ./ p.T1))
            end
            ADC_idx += 1
        end

        Bz_1, Bz_2 = Bz_2, Bz_1
    end

    #Final Spin-State
    M.xy .= M.xy .* exp.(-t_seq ./ p.T2) .* (cos.(ϕ) .+ im * sin.(ϕ))
    M.z .= M.z .* exp.(-t_seq ./ p.T1) .+ p.ρ .* (1 .- exp.(-t_seq ./ p.T1))

    return nothing
end
