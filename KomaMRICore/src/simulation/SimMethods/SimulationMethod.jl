"""
    Bloch()

Default Bloch simulation using the hard-pulse approximation: the magnetic field
`B` [T] is treated as piecewise constant over each integration step `Δt` [s].
"""
struct Bloch <: SimulationMethod end
export Bloch

include("BlochMagnus/BlochMagnus.jl")

include("Magnetization.jl")

@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl

function sim_output_dim(
    obj::Phantom, seq::Sequence, sys::Scanner, sim_method::SimulationMethod
)
    return (sum(seq.ADC.N), 1) #Nt x Ncoils, This should consider the coil info from sys
end

function split_sig_per_thread(sig, i, p, sim_method::SimulationMethod)
    return @view sig[:, :, i]
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom, sim_method::SimulationMethod)
    Mxy = complex.(zero.(obj.ρ))
    Mz = copy(obj.ρ)
    Xt = Mag(Mxy, Mz)
    return Xt, obj
end

"""Stores pre-allocated arrays for use in run_spin_precession! and run_spin_excitation!"""
abstract type PreallocResult end

"""Default preallocation struct, stores nothing."""
struct DefaultPrealloc <: PreallocResult end

Base.view(p::PreallocResult, i::UnitRange) = p

"""Default preallocation function."""
prealloc(sim_method::SimulationMethod, backend::KA.Backend, obj::Phantom, M::Mag, max_block_length::Integer, groupsize) = DefaultPrealloc()

include("BlochSimple/BlochSimple.jl")
include("Bloch/cpu/BlochCPU.jl")
include("BlochMagnus/cpu/BlochMagnusCPU.jl")
include("Bloch/gpu/BlochGPU.jl")
include("BlochMagnus/gpu/MagnusMidKernel.jl")
include("BlochMagnus/gpu/MagnusQuadKernel.jl")
include("BlochMagnus/gpu/MagnusGLKernel.jl")
include("BlochMagnus/gpu/MagnusBGLKernel.jl")
include("BlochDict/BlochDict.jl")
