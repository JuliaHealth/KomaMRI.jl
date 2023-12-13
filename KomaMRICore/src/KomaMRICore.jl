module KomaMRICore

using KomaMRIBase
#General
import Base.*, Base.abs
using ThreadsX
#Printing
using ProgressMeter
#Simulation
using CUDA

#Rawdata
include("rawdata/ISMRMRD.jl")
#Datatypes
include("datatypes/Spinor.jl")
include("datatypes/Magnetization.jl")
include("other/DiffusionModel.jl")
#Simulator
include("simulation/GPUFunctions.jl")
include("simulation/SimulatorCore.jl")

#ISMRMRD
export signal_to_raw_data
#Simulator
export simulate, simulate_slice_profile

#Package version, KomaMRICore.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
