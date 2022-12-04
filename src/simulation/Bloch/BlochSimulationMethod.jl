struct Bloch <: SimulationMethod end
include("Magnetization.jl") #Defines Mag <: SpinStateRepresentation
@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom{T}, sim_method::Bloch) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt
end

"""
    run_spin_precession(obj, seq, Xt, sig)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∫ ρ(x,t) exp(- t/T2(x,t) ) exp(- 𝒊 ϕ(x,t)) dx. It performs the simulation in free
precession.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `M0`: (`::Vector{Mag}`) initial state of the Mag vector (actually, it's a part of the
    complete Mag vector)

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
NVTX.@range function run_spin_precession(p::Phantom{T}, s::DiscreteSequence{T}, M::Mag{T}) where {T<:Real}
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, s.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, s.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, s.t')
    #Effective field
    Bz = xt .* s.Gx' .+ yt .* s.Gy' .+ zt .* s.Gz'
    #Rotation
    if is_ADC_on(s)
        ϕ = T(2π * γ) .* cumtrapz(s.Δt', Bz)
    else
        ϕ = T(2π * γ) .* trapz(s.Δt', Bz)
    end
    #Mxy preccesion and relaxation and Mz relaxation
    tp = cumsum(s.Δt) # t' = t - t0
    dur = sum(s.Δt) #Total length, used for signal relaxation
    Mxy = M.xy .* exp.(1im .* (ϕ .+ p.Δw .* tp') .- tp' ./ p.T2) #This assumes Δw and T2 are constant
    M.xy .= Mxy[:, end]
    M.z .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    #Acquired signal
    sig = sum(Mxy[:, findall(s.ADC)]; dims=1) #<--- TODO: add coil sensitivities
    return sig, M
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
NVTX.@range function run_spin_excitation(p::Phantom{T}, seq::DiscreteSequence{T}, M::Mag{T}) where {T<:Real}
    #Simulation
    for s ∈ seq
        #Motion
        xt = p.x .+ p.ux(p.x, p.y, p.z, s.t)
        yt = p.y .+ p.uy(p.x, p.y, p.z, s.t)
        zt = p.z .+ p.uz(p.x, p.y, p.z, s.t)
        #Effective field
        ΔB0 = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
        Bz = (s.Gx .* xt .+ s.Gy .* yt .+ s.Gz .* zt) .+ ΔB0 #<-- TODO: This line is very slow, FIX!?
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here,  this is just Forward Euler
        M = Q(φ, s.B1 ./ B, Bz ./ B) * M
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z .= M.z .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    return M
end
