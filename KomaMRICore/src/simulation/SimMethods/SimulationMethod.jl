#Bloch is the default simulation method
struct Bloch <: SimulationMethod end

export Bloch
include("Magnetization.jl")

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

abstract type PreallocResult{T<:Real} end

"""Default preallocation struct, stores nothing."""
struct DefaultPreAlloc{T} <: PreallocResult{T} end

Base.view(p::DefaultPreAlloc, i::UnitRange) = p

"""Default preallocation function."""
prealloc(sim_method::SimulationMethod, backend::KA.Backend, obj::Phantom{T}, M::Mag{T}) where {T<:Real} = DefaultPreAlloc{T}()

include("../../datatypes/VectorSU2.jl")
include("KernelFunctions.jl")
include("BlochSimple/BlochSimple.jl")
include("Bloch/BlochCPU.jl")
include("BlochDict/BlochDict.jl")