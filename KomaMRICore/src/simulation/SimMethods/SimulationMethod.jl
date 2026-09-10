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
    return (sum(seq.ADC.N), get_n_coils(sys.receiver))
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

# Coupled motions must compute all three displacements before modifying positions.
struct MotionCoordinates{A}
    positions::NTuple{3,A}
    displacements::NTuple{3,A}
end

function Base.view(
    coordinates::MotionCoordinates{<:AbstractVector}, i,
)
    return MotionCoordinates(
        map(a -> view(a, i), coordinates.positions),
        map(a -> view(a, i), coordinates.displacements),
    )
end

view_motion_coordinates(motion::NoMotion, _) = motion
view_motion_coordinates(coordinates::MotionCoordinates, i) = view(coordinates, i)
motion_enabled(::NoMotion) = Val(false)
motion_enabled(::MotionCoordinates) = Val(true)

prealloc_motion_coordinates(motion::NoMotion, _, _, _) = motion
function prealloc_motion_coordinates(
    ::Union{Motion,MotionList}, ::KA.CPU, obj, _,
)
    buffers() = ntuple(_ -> similar(obj.x), 3)
    return MotionCoordinates(buffers(), buffers())
end
function prealloc_motion_coordinates(
    ::Union{Motion,MotionList}, backend::KA.GPU, obj, max_block_length,
)
    T = eltype(obj.x)
    buffers() = ntuple(_ -> KA.zeros(backend, T, length(obj), max_block_length), 3)
    return MotionCoordinates(buffers(), buffers())
end

spin_coordinates!(::NoMotion, _, x, y, z, _) = (x, y, z)
function spin_coordinates!(
    coordinates::MotionCoordinates{<:AbstractVector}, motion, x, y, z, t,
)
    return KomaMRIBase.get_spin_coords!(
        coordinates.positions, coordinates.displacements, motion, x, y, z, t,
    )
end
function spin_coordinates!(
    coordinates::MotionCoordinates{<:AbstractMatrix}, motion, x, y, z, t,
)
    columns = 1:length(t)
    positions = map(a -> view(a, :, columns), coordinates.positions)
    displacements = map(a -> view(a, :, columns), coordinates.displacements)
    return KomaMRIBase.get_spin_coords!(
        positions, displacements, motion, x, y, z, t,
    )
end

"""Default preallocation function."""
prealloc(
    ::SimulationMethod, ::KA.Backend, obj, _M,
    _max_block_length, _max_adc_samples, _groupsize, _sys,
) = DefaultPrealloc()

include("BlochSimple/BlochSimple.jl")
include("Bloch/cpu/BlochCPU.jl")
include("BlochMagnus/cpu/BlochMagnusCPU.jl")
include("Bloch/gpu/BlochGPU.jl")
include("BlochMagnus/gpu/MagnusMidKernel.jl")
include("BlochMagnus/gpu/MagnusQuadKernel.jl")
include("BlochMagnus/gpu/MagnusGLKernel.jl")
include("BlochMagnus/gpu/MagnusBGLKernel.jl")
include("BlochDict/BlochDict.jl")
