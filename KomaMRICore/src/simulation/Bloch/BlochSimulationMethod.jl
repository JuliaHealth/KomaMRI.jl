struct Bloch <: SimulationMethod end

export Bloch

include("Magnetization.jl") #Defines Mag <: SpinStateRepresentation
@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl

function sim_output_dim(
    obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::SimulationMethod
) where {T<:Real}
    return (sum(seq.ADC.N), 1) #Nt x Ncoils, This should consider the coil info from sys
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(
    obj::Phantom{T}, sim_method::SimulationMethod
) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
end

"""Preallocated arrays for use in run_spin_precession."""
struct BlochPrealloc{T} <: PreallocResult{T}
    Bz_old::AbstractVector{T}
    Bz_new::AbstractVector{T}
    ϕ::AbstractVector{T}
    Mxy::AbstractVector{Complex{T}}
end

Base.view(p::BlochPrealloc, i::UnitRange) = begin
    @views BlochPrealloc(
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Mxy[i]
    )
end

"""Default preallocation function. Returns arrays for use in run_spin_precession."""
function prealloc(sim_method::SimulationMethod, obj::Phantom{T}, M::Mag{T}) where {T<:Real}
    BlochPrealloc(
        similar(obj.x),
        similar(obj.x),
        similar(obj.x),
        similar(M.xy)
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
    sim_method::SimulationMethod,
    backend::KA.Backend,
    prealloc::BlochPrealloc
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
function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::SimulationMethod,
) where {T<:Real}
    #Simulation
    for s in seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, s.t)
        #Effective field
        ΔBz = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        Bz = (s.Gx .* x .+ s.Gy .* y .+ s.Gz .* z) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        mul!(Q(φ, s.B1 ./ B, Bz ./ B), M)
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z .= M.z .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    #Acquired signal
    #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end
