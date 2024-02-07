using Pkg
Pkg.activate(".")

using KomaMRI  

KomaMRICore.CUDA.allowscalar() = false;

phantom = brain_phantom2D()[1:100]

phantom.ux = (x,y,z,t) -> sin.(2π*t)

sys = Scanner()

seq = PulseDesigner.EPI_example(; sys)

## Simulation
simulate(phantom, seq, sys)

simulate(phantom, seq, sys)