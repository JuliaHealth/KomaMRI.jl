module KomaMRICore

#IMPORT PACKAGES
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty
#General
using InteractiveUtils      # for using subtypes()
using Reexport, ThreadsX
#Printing
using Scanf, ProgressMeter
#Datatypes
using Parameters
#Simulation
using Interpolations
using CUDA
#Reconstruction
using MRIBase
@reexport using MRIBase: EncodingCounters, Profile, RawAcquisitionData, AcquisitionData, AcquisitionHeader
using MAT   # For loading example phantoms

global γ = 42.5774688e6; #Hz/T gyromagnetic constant for H1, JEMRIS uses 42.5756 MHz/T

#Hardware
include("datatypes/Scanner.jl")
#Sequence
include("datatypes/sequence/Grad.jl")
include("datatypes/sequence/RF.jl")
include("datatypes/sequence/ADC.jl")
include("simulation/KeyValuesCalculation.jl")
include("datatypes/Sequence.jl")
include("datatypes/sequence/Delay.jl")
#Phantom
include("datatypes/Phantom.jl")
#Rawdata
include("rawdata/ISMRMRD.jl")
#Simulator
include("datatypes/simulation/DiscreteSequence.jl")
include("datatypes/simulation/Spinor.jl")
include("simulation/TimeStepCalculation.jl")
include("simulation/other/DiffusionModel.jl")
include("simulation/GPUFunctions.jl")
# include("simulation/other/OffResonanceModel.jl")
include("simulation/TrapezoidalIntegration.jl")
include("simulation/SimulatorCore.jl")

#Main
export γ #gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export Mag, dur
#Phantom
export brain_phantom2D, brain_phantom3D
#Spinors
export Spinor, Rx, Ry, Rz, Q, Un
#Secondary
export get_kspace, rotx, roty, rotz
#ISMRMRD
export signal_to_raw_data
#Simulator
export simulate, simulate_slice_profile

#Additionals
export get_flip_angles, is_RF_on, is_GR_on, is_ADC_on
#Sequence related
export get_M0, get_M1, get_M2, get_kspace

#PulseDesigner submodule
include("sequences/PulseDesigner.jl")
export PulseDesigner

#Package version, KomaMRICore.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
