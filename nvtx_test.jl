using Pkg
Pkg.activate(".")

using KomaMRI  

KomaMRICore.CUDA.allowscalar() = false;

phantom = brain_phantom2D()[1:100]

phantom.ux = (x,y,z,t)->sin.(x.*t)
phantom.uy = (x,y,z,t)->sin.(t)
phantom.uz = (x,y,z,t)->0

sys = Scanner()

seq = PulseDesigner.EPI_example(; sys)

## Simulation
simulate(phantom, seq, sys)

simulate(phantom, seq, sys)