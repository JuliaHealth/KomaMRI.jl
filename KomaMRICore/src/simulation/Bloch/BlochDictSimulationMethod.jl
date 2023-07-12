Base.@kwdef struct BlochDict <: SimulationMethod 
    save_Mz::Bool=false
end

export BlochDict
Base.show(io::IO, s::BlochDict) = begin
	print(io, "BlochDict(save_Mz=$(s.save_Mz))")
end


output_Ndim(sim_method::BlochDict) = 3

function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochDict) where {T<:Real}
    out_state_dim = sim_method.save_Mz ? 2 : 1
    return (sum(seq.ADC.N), length(obj), out_state_dim)
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom{T}, sim_method::BlochDict) where {T<:Real}
    return initialize_spins_state(obj, Bloch())
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
function run_spin_precession!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
    M::Mag{T}, sim_method::BlochDict) where {T<:Real}
    #Simulation
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    #Effective field
    Bz = xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz' .+ p.Δw / T(2π * γ)
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π * γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π * γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy precession and relaxation, and Mz relaxation
    tp = cumsum(seq.Δt) # t' = t - t0
    dur = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy = [M.xy M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2)] #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    
    #Acquired signal
    sig[:,:,1] .= transpose(Mxy[:, findall(seq.ADC)])

    if sim_method.save_Mz
        Mz = [M.z M.z .* exp.(-tp' ./ p.T1) .+ p.ρ .* (1 .- exp.(-tp' ./ p.T1))] #Calculate intermediate points
        sig[:,:,2] .= transpose(Mz[:, findall(seq.ADC)]) #Save state to signal
        M.z .= Mz[:, end]
    else
        M.z .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1)) #Jump to the last point
    end
    return nothing
end

"""
    M0 = run_spin_excitation(obj, seq, M0)

It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (actually, it's
    a part of the complete Mag vector and it's a part of the initial state for the next
    precession simulation step)
"""
function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochDict) where {T<:Real}
    run_spin_excitation!(p, seq, sig, M, Bloch()) #The same as Bloch
    return nothing
end
