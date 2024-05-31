module KomaMRICore

# General
import Base.*, Base.abs
using Reexport
using ThreadsX
# Printing
using ProgressMeter
# Simulation
using CUDA

# KomaMRIBase
@reexport using KomaMRIBase
@reexport import KomaMRIBase.get_spin_coords # This should not be necessary, but it is

# Rawdata
include("rawdata/ISMRMRD.jl")
# Datatypes
include("datatypes/Spinor.jl")
include("other/DiffusionModel.jl")
# Simulator
include("simulation/GPUArbitraryMotion.jl")
include("simulation/GPUFunctions.jl")
include("simulation/SimulatorCore.jl")

# ISMRMRD
export signal_to_raw_data
# Simulator
export Mag
export simulate, simulate_slice_profile
# Spin coordinates
export get_spin_coords
# Spinors
export Spinor, Rx, Ry, Rz, Q, Un

#Package version, KomaMRICore.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
