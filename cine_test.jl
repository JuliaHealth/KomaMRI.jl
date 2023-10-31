using KomaMRI

include("CardiacCine.jl")

hr = 60         # [bpm]
N = 32          # image size = N x N
N_phases = 5    # Number of cardiac phases
 
# ------------ CHOOSE ONE OF THIS 2 PHANTOMS: -----------
## Left Ventricle Phantom 
path = "examples/2.phantoms/ring_motion.phantom"
FOV = 0.15

## Flow Phantom
path = "examples/2.phantoms/foo_flow.phantom"
FOV = 0.03    
## ------------------------------------------------------

phantom = read_phantom(path)

sys = Scanner()

frames = cardiac_cine(FOV,hr,N_phases,N,phantom,sys);

plot_cine(frames,N_phases)