module KomaMRI

#IMPORT PACKAGES
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size,
       Base.copy, Base.Threads.@spawn, Base.Threads.@threads,
       Base.angle, Base.abs, Base.getproperty, Base.one, Base.zero
#General
using Pkg, Random, Reexport, ThreadsX
#Printing
using Scanf, ProgressMeter
#Datatypes
using Parameters
#Simulation
using CUDA, Interpolations, Hwloc
#Reconstruction
using MRIReco, MRIFiles
@reexport using MRIReco: RawAcquisitionData, AcquisitionData, reconstruction
@reexport using MRIFiles: ISMRMRDFile
#IO
using FileIO, HDF5, MAT, JLD2
#GUI
using Blink, Interact, PlotlyJS, AssetRegistry
@reexport using PlotlyJS: savefig

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
include("sequences/PulseDesigner.jl")
include("io/Pulseq.jl")
#Phantom
include("datatypes/Phantom.jl")
include("io/JEMRIS.jl")
include("io/MRiLab.jl")
#Reconstruction
include("reconstruction/Recon.jl")
#Simulator
include("datatypes/simulation/Magnetization.jl")
include("simulation/TimeStepCalculation.jl")
include("simulation/other/DiffusionModel.jl")
# include("simulation/other/OffResonanceModel.jl")
include("simulation/TrapezoidalIntegration.jl")
include("simulation/SimulatorCore.jl")
include("io/ISMRMRD.jl")
#UI
include("ui/DisplayFunctions.jl")

#Main
export γ #gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export Mag, dur
#Pulseq
export read_seq
#ISMRMRD
export signal_to_raw_data
# Phantom
export brain_phantom2D, brain_phantom3D, read_phantom_jemris, read_phantom_MRiLab
#RF-related
export Spinor, Rx, Ry, Rz, Q, Un
#Secondary
export PulseDesigner, get_kspace, rotx, roty, rotz
# Display
export plot_seq, plot_grads_moments, plot_kspace, plot_phantom_map, plot_signal, plot_M0, plot_image
# Simulator
export simulate, simulate_slice_profile
#GUI
!Blink.AtomShell.isinstalled() && Blink.AtomShell.install()
include("KomaUI.jl")

export KomaUI

end
