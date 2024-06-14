module KomaMRICore

# General
import Base.*, Base.abs
import KernelAbstractions as KA
using Reexport
using ThreadsX
# Printing
using ProgressMeter

# KomaMRIBase
@reexport using KomaMRIBase

# Rawdata
include("rawdata/ISMRMRD.jl")
# Datatypes
include("datatypes/Spinor.jl")
include("other/DiffusionModel.jl")
# Simulator
include("simulation/Functors.jl")
include("simulation/GPUFunctions.jl")
include("simulation/SimulatorCore.jl")

# ISMRMRD
export signal_to_raw_data
# Simulator
export Mag
export simulate, simulate_slice_profile
# Spinors
export Spinor, Rx, Ry, Rz, Q, Un

#Package version, KomaMRICore.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
