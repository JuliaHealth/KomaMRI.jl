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
    run_spin_precession!(obj, seqd, sig, M, sim_method)

Conduct the simulation within the precession regime using the `BlochDict` simulation method.
It computes the magnetization in the xy plane and in z axis. The raw signal `sig` and the
magnetization state `M` are updated in-place, representing the result of the simulation.

# Arguments
- `obj`: (`::Phantom{T:<Real}`) Phantom struct
- `seqd`: (`::DiscreteSequence{T:<Real}`) DiscreteSequence struct
- `sig`: (`::AbstractArray{Complex{T:<Real}}`) raw signal
- `M`: (`::Mag{T:<Real}`) magnetization state
- `sim_method`: (`::BlochDict`) utilized for dispatching the `BlochDict` simulation method
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
    run_spin_excitation!(obj, seqd, sig, M, sim_method)

Conduct the simulation within the excitation regime using the `BlochDict` simulation method.
It computes the same as the `Bloch` simulation method. The raw signal `sig` and the
magnetization state `M` are updated in-place, representing the result of the simulation.

# Arguments
- `obj`: (`::Phantom{T:<Real}`) Phantom struct
- `seqd`: (`::DiscreteSequence{T:<Real}`) DiscreteSequence struct
- `sig`: (`::AbstractArray{Complex{T:<Real}}`) raw signal
- `M`: (`::Mag{T:<Real}`) magnetization state
- `sim_method`: (`::BlochDict`) utilized for dispatching the `BlochDict` simulation method
"""
function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochDict) where {T<:Real}
    run_spin_excitation!(p, seq, sig, M, Bloch()) #The same as Bloch
    return nothing
end
