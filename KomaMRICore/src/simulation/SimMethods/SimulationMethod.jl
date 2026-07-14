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
    return (sum(seq.ADC.N), get_n_coils(sys.receiver)) #Nt x Ncoils, This should consider the coil info from sys
end

function split_sig_per_thread(sig, i, p, sim_method::SimulationMethod)
    return @view sig[:, :, i]
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(
    obj::Phantom{T}, sim_method::SimulationMethod
) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(Complex{T}, Nspins)
    Mz = copy(obj.ρ)
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
end

"""Stores pre-allocated arrays for use in run_spin_precession! and run_spin_excitation!"""
abstract type PreallocResult{T<:Real} end

"""Default preallocation struct, stores nothing."""
struct DefaultPrealloc{T} <: PreallocResult{T} end

Base.view(p::PreallocResult, i::UnitRange) = p

"""Default preallocation function."""
prealloc(sim_method::SimulationMethod, backend::KA.Backend, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize, sys::Scanner) where {T<:Real} = DefaultPrealloc{T}()

include("BlochSimple/BlochSimple.jl")
include("Bloch/cpu/BlochCPU.jl")
include("BlochMagnus/cpu/BlochMagnusCPU.jl")
include("Bloch/gpu/BlochGPU.jl")
include("BlochMagnus/gpu/MagnusMidKernel.jl")
include("BlochMagnus/gpu/MagnusQuadKernel.jl")
include("BlochMagnus/gpu/MagnusGLKernel.jl")
include("BlochMagnus/gpu/MagnusBGLKernel.jl")
include("BlochDict/BlochDict.jl")
