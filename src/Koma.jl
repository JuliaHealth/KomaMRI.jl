module Koma

#IMPORT PACKAGES
using LinearAlgebra: Matrix
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, 
       Base.copy, Base.Threads.@spawn, Base.Threads.@threads,
       Base.angle, Base.abs, Base.getproperty, Base.one
using Random, LinearAlgebra, FFTW, Images, Printf, MAT, PlotlyJS, 
       ProgressMeter, CUDA, Parameters, Interpolations, FileIO

global γ = 42.5774688e6; #Hz/T gyromagnetic constant for H1

#Hardware
include("datatypes/Scanner.jl")

#Sequence
include("datatypes/sequence/Grad.jl")
include("datatypes/sequence/RF.jl")
include("datatypes/sequence/ADC.jl")
include("datatypes/Sequence.jl")
include("datatypes/sequence/Delay.jl")
include("sequences/PulseDesigner.jl")
include("io/Pulseq.jl")

#Phantom
include("datatypes/Phantom.jl")
#Reconstruction
include("reconstruction/Recon.jl")
#Simulator
include("datatypes/simulation/Magnetization.jl")
# include("simulation/DiffusionModel.jl")
# include("simulation/OffResonanceModel.jl")
include("simulation/Simulator.jl")
#UI
include("ui/Display.jl")

#Main 
export γ #gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export Grad, RF, ADC, Delay
export Mag, dur
#RF-related
export Spinor, Rx, Ry, Rz, Q, Un
#Secondary
export PulseDesigner, get_designed_kspace, rotx, roty, rotz
# Display
export plot_seq, plot_grads_moments, plot_ksapce_trajectory
# Simulator
export simulate, run_sim_time_iter

#GUI
using Blink, Interact, AssetRegistry, JLD2
!Blink.AtomShell.isinstalled() && Blink.AtomShell.install()
include("KomaUI.jl")

export KomaUI

end
