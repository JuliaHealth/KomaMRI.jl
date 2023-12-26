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
include("other/DiffusionModel.jl")
#Simulator
include("simulation/GPUFunctions.jl")
include("simulation/SimulatorCore.jl")

#Main
export γ #gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export Mag, dur
#Phantom
export brain_phantom2D, brain_phantom3D
#Secondary
export get_kspace, rotx, roty, rotz
#ISMRMRD
export signal_to_raw_data
#Simulator
export simulate, simulate_slice_profile

#Package version, KomaMRICore.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
