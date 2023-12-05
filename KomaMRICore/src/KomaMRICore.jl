module KomaMRICore

using KomaMRIBase
#General
using ThreadsX
#Printing
using Scanf, ProgressMeter
#Simulation
using CUDA

#Rawdata
include("rawdata/ISMRMRD.jl")
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
