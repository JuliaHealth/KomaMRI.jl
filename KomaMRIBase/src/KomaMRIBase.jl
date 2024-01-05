module KomaMRIBase

#IMPORT PACKAGES
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty
#General
using Reexport
#Datatypes
using Parameters
#Simulation
using Interpolations
#Reconstruction
using MRIBase
@reexport using MRIBase: Profile, RawAcquisitionData, AcquisitionData, AcquisitionHeader, EncodingCounters, Limit
using MAT   # For loading example phantoms

global γ = 42.5774688e6; # Hz/T gyromagnetic constant for H1, JEMRIS uses 42.5756 MHz/T

# Hardware
include("datatypes/Scanner.jl")
# Sequence
include("datatypes/sequence/Grad.jl")
include("datatypes/sequence/RF.jl")
include("datatypes/sequence/ADC.jl")
include("timing/KeyValuesCalculation.jl")
include("datatypes/Sequence.jl")
include("datatypes/sequence/Delay.jl")
# Phantom
include("datatypes/Phantom.jl")
# Simulator
include("datatypes/simulation/DiscreteSequence.jl")
include("timing/TimeStepCalculation.jl")
include("timing/TrapezoidalIntegration.jl")

# Main
export γ    # gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export dur, get_block_start_times, get_samples
export DiscreteSequence
export discretize, get_adc_phase_compensation, get_adc_sampling_times
# This are also used for simulation
export kfoldperm, trapz, cumtrapz
# Phantom
export brain_phantom2D, brain_phantom3D, pelvis_phantom2D
# Spinors
export Spinor, Rx, Ry, Rz, Q, Un
# Secondary
export get_kspace, rotx, roty, rotz
# Additionals
export get_flip_angles, is_RF_on, is_GR_on, is_ADC_on
# Sequence related
export get_M0, get_M1, get_M2, get_kspace

# PulseDesigner submodule
include("sequences/PulseDesigner.jl")
export PulseDesigner

#Package version, KomaMRIBase.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaMRIBase
