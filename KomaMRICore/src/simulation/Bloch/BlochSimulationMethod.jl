struct Bloch <: SimulationMethod end

export Bloch

include("Magnetization.jl") #Defines Mag <: SpinStateRepresentation
@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl

output_Ndim(sim_method::Bloch) = 2 #time-points x coils

function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::Bloch) where {T<:Real}
    return (sum(seq.ADC.N), 1) #Nt x Ncoils, This should consider the coil info from sys
end

"""
Magnetization initialization for Bloch simulation method.
"""
function initialize_spins_state(obj::Phantom{T}, sim_method::Bloch) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
end

"""
    run_spin_precession!(obj, seqd, sig, M, sim_method)

Conduct the simulation within the precession regime using the `Bloch` simulation method. It
calculates S(t) = ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). The raw signal `sig` and
the magnetization state `M` are updated in-place, representing the result of the simulation.

# Arguments
- `obj`: (`::Phantom{T:<Real}`) Phantom struct
- `seqd`: (`::DiscreteSequence{T:<Real}`) DiscreteSequence struct
- `sig`: (`::AbstractArray{Complex{T:<Real}}`) raw signal
- `M`: (`::Mag{T:<Real}`) magnetization state
- `sim_method`: (`::Bloch`) utilized for dispatching the `Block` simulation method
"""
function run_spin_precession!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
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
    M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    #Acquired signal
    sig .= transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1)) #<--- TODO: add coil sensitivities
    return nothing
end

"""
    run_spin_excitation!(obj, seqd, sig, M, sim_method)

Conduct the simulation within the excitation regime using the `Bloch` simulation method.
This regime gives rise to a rotation of magnetization with an angle determined by the
effective magnetic field (including B1, gradients, and off-resonance) and with respect to a
rotation axis. The raw signal `sig` and the magnetization state `M` are updated in-place,
representing the result of the simulation.

# Arguments
- `obj`: (`::Phantom{T:<Real}`) Phantom struct
- `seqd`: (`::DiscreteSequence{T:<Real}`) DiscreteSequence struct
- `sig`: (`::AbstractArray{Complex{T:<Real}}`) raw signal
- `M`: (`::Mag{T:<Real}`) magnetization state
- `sim_method`: (`::Bloch`) utilized for dispatching the `Block` simulation method
"""
function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    #Simulation
    for s ∈ seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        xt = p.x .+ p.ux(p.x, p.y, p.z, s.t)
        yt = p.y .+ p.uy(p.x, p.y, p.z, s.t)
        zt = p.z .+ p.uz(p.x, p.y, p.z, s.t)
        #Effective field
        ΔBz = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
        Bz = (s.Gx .* xt .+ s.Gy .* yt .+ s.Gz .* zt) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        mul!( Q(φ, s.B1 ./ B, Bz ./ B), M )
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z  .= M.z  .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    #Acquired signal
    #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end
