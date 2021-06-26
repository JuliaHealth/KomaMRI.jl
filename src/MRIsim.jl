module MRIsim

#IMPORT PACKAGES
using LinearAlgebra: Matrix
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, 
       Base.copy, Base.Threads.@spawn, Base.Threads.@threads,
       Base.angle, Base.abs, Base.getproperty
using Random, LinearAlgebra, FFTW, Images, Printf, MAT, PlotlyJS, 
       ProgressMeter, CUDA, Parameters

global γ = 42.5e6; #Hz/T gyromagnetic constant for H1

#CORE
include("Grad.jl")
include("RF.jl")
include("Sequence.jl")
include("Phantom.jl")
include("Magnetization.jl")
include("Simulator.jl")
include("Recon.jl")
include("Display.jl")

#UNDER DEVELOPMENT
include("DiffusionModel.jl")

#Main 
export Grad, delay, RF, Sequence, Phantom, Mag
#RF-related
export Rx, Ry, Rz, Q, Un
#Secondary
export PulseDesigner

#GUI
using Blink, Interact, AssetRegistry, JLD2, FileIO
!Blink.AtomShell.isinstalled() && Blink.AtomShell.install()
include("SpinLab.jl")

export SpinLab

end
