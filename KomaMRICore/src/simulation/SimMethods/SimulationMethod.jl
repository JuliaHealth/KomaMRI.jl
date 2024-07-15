#Bloch is the default simulation method
struct Bloch <: SimulationMethod end

export Bloch

include("Magnetization.jl")

@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl

function sim_output_dim(
    obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::SimulationMethod
) where {T<:Real}
    # Determine the number of coils
    n_coils = size(obj.coil_sens, 2)
    return (sum(seq.ADC.N), n_coils) # Nt x Ncoils
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

"""Stores pre-allocated arrays for use in run_spin_precession! and run_spin_excitation!"""
abstract type PreallocResult{T<:Real} end

"""Stores information precalculated before the simulation objects are moved to the GPU."""
abstract type PrecalcResult{T<:Real} end

"""Default preallocation struct, stores nothing."""
struct DefaultPrealloc{T} <: PreallocResult{T} end

Base.view(p::PreallocResult, i::UnitRange) = p
prealloc_block(p::PreallocResult, i::Integer) = p

"""Default precalculation struct, stores nothing."""
struct DefaultPrecalc{T} <: PrecalcResult{T} end

"""Default preallocation function."""
prealloc(sim_method::SimulationMethod, backend::KA.Backend, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, precalc) where {T<:Real} = DefaultPrealloc{T}()

"""Default precalc function."""
precalculate(sim_method::SimulationMethod, backend::KA.Backend, seq::DiscreteSequence{T}, parts::Vector{UnitRange{S}}, excitation_bool::Vector{Bool}) where {T<:Real,S<:Integer} = DefaultPrecalc{T}()

include("BlochSimple/BlochSimple.jl")
include("Bloch/BlochCPU.jl")
include("Bloch/BlochGPU.jl")
include("BlochDict/BlochDict.jl")