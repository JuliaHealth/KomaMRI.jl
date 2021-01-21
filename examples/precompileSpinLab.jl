# This is an example code used to precompile the SpinLab GUI
# Using a precompiled version of MRIsim should accelerate startup times for each library
# In the future, this will enable to create an Standalone-App (exceutable)
using MRIsim
SpinLab()

##
phantom = MRIsim.brain_phantom2D(;axis="coronal")
EPI,_,_,_ = MRIsim.EPI_base(40/100, 100, 4e-6, 60e-3)
TE = 25e-3 
d = MRIsim.delay(TE-MRIsim.dur(EPI)/2)
DELAY = Sequence([d;d])
seq = DELAY + EPI
println("EPI successfully loaded! (TE = $(TE*1e3) ms)")
scanner = []
signal = 0
kdata = [0.0im 0.; 0. 0.]
MRIsim.print_gpus_info()

Δt = 4e-6 #<- simulate param
t = collect(Δt:Δt:MRIsim.dur(seq))
Nphant, Nt = prod(size(phantom)), length(t)
N_parts = floor(Int, Nphant*Nt/2.7e6*1/3) #heuristic
println("Dividing simulation in Nblocks=$N_parts")
S = @time MRIsim.run_sim2D_times_iter(phantom,seq,t;N_parts) #run_sim2D_times_iter run_sim2D_spin
signal = S[MRIsim.get_DAC_on(seq,t)]/prod(size(phantom)) #Acquired data
S = nothing
Nx = Ny = 100 #hardcoded by now
kdata = reshape(signal,(Nx,Ny)) #Turning into kspace image
kdata[:,2:2:Ny,:] = kdata[Nx:-1:1,2:2:Ny] #Flip in freq-dir for the EPI
kdata = convert(Array{Complex{Float64},2},kdata);
##
using CUDA
a = cu([1,2])
b = a .* a

@info "GPU related operations do not work yet, wait until Julia 1.6"

##
using PackageCompiler
PackageCompiler.create_sysimage(:MRIsim;
sysimage_path="C:/Users/ccp/Documents/MRIsim.jl/src/examples/sysimageSpinLab.dll",
precompile_execution_file="C:/Users/ccp/Documents/MRIsim.jl/src/examples/precompileSpinLab.jl")
# julia --sysimage C:/Users/ccp/Documents/MRIsim.jl/src/examples/sysimageSpinLab.dll --project=@.

##
# using PackageCompiler
# PackageCompiler.create_sysimage([:Plots, :JuMP, :Ipopt, :Blink, :PlotlyJS], replace_default=true)