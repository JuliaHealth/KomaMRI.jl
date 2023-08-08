module KomaMRICore

#IMPORT PACKAGES
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty
#General
using Reexport, ThreadsX
#Printing
using Scanf, ProgressMeter
#Datatypes
using Parameters
#Simulation
using Interpolations
#Reconstruction
using MRIBase, MRIFiles
@reexport using MRIBase: Profile, RawAcquisitionData, AcquisitionData, AcquisitionHeader
@reexport using MRIFiles: ISMRMRDFile
#IO
using FileIO, HDF5, MAT, JLD2

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
include("io/Pulseq.jl")
#Phantom
include("datatypes/Phantom.jl")
include("io/JEMRIS.jl")
include("io/MRiLab.jl")
#Simulator
include("datatypes/simulation/DiscreteSequence.jl")
include("datatypes/simulation/Spinor.jl")
include("simulation/TimeStepCalculation.jl")
include("simulation/other/DiffusionModel.jl")
include("simulation/GPUFunctions.jl")
# include("simulation/other/OffResonanceModel.jl")
include("simulation/TrapezoidalIntegration.jl")
include("simulation/SimulatorCore.jl")
include("io/ISMRMRD.jl")

#Main
export γ #gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export Mag, dur
#Pulseq
export read_seq
#ISMRMRD
export signal_to_raw_data
#Phantom
export brain_phantom2D, brain_phantom3D, read_phantom_jemris, read_phantom_MRiLab
#Spinors
export Spinor, Rx, Ry, Rz, Q, Un
#Secondary
export get_kspace, rotx, roty, rotz
#Simulator
export simulate, simulate_slice_profile

#Additionals
export get_flip_angles, is_RF_on, is_GR_on, is_ADC_on
#Sequence related
export get_M0, get_M1, get_M2, get_kspace

#PulseDesigner submodule
include("sequences/PulseDesigner.jl")
export PulseDesigner

using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

using PackageExtensionCompat
function __init__()
    @require_extensions
end
#PrecompileTools
# using PrecompileTools
# @setup_workload begin
#     obj = brain_phantom2D();
#     sys = Scanner();
#     seq = read_seq(joinpath(dirname(pathof(KomaMRICore)), 
#         "../../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq"));
#     simParams = KomaMRICore.default_sim_params();
#     simParams["gpu"] = true;
#     @compile_workload begin
#         raw = simulate(obj, seq, sys; simParams);
#         #abs.(raw.profiles[1].data)
#     end
# end

end
