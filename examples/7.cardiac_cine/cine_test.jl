using KomaMRI

include("CardiacCine.jl")

hr = 60         # [bpm]
N = 64          # image size = N x N
N_phases = 6    # Number of cardiac phases
 
# ------------ CHOOSE ONE OF THIS PHANTOMS: -----------
## Left Ventricle Phantom 
path = "examples/2.phantoms/ring_motion.phantom"

## Flow Phantom
path = "examples/2.phantoms/foo_flow.phantom"   

## Flow vertical
path = "examples/2.phantoms/flow_artery_vertical.phantom" 

## Flow horizontal
path = "examples/2.phantoms/flow_artery_horizontal.phantom" 

##
phantom = read_phantom(path)

## Motion heart
phantom = read_phantom_MAT("/datos/work/phantomXCAT_1mm_5cardPhases/"; ss=1, Δx=1)
phantom = phantom[abs.(phantom.z) .< 0.005]
phantom = phantom[abs.(phantom.x) .<= 0.06]
phantom = phantom[abs.(phantom.y) .<= 0.06]
phantom.motion.Δz = zeros(length(phantom),phantom.motion.K - 1)
## ------------------------------------------------------

FOV = 2*maximum([maximum(abs.(phantom.x)),
                 maximum(abs.(phantom.y)),
                 maximum(abs.(phantom.z))])

sys = Scanner()

FOV = 0.15
TR=80e-3
## Simulate
frames = cardiac_cine(FOV,hr,N_phases,N,phantom,sys;
                      Δf=0,
                      flip_angle=45, 
                      TR=TR,
                      dummy_cycles = 0,
                      tagging=false);

## Post-processing
post = []
post2 = []
for frame in frames
    aux = copy(frame)
    perc = KomaMRICore.percentile(frame[:],96)
    aux[frame.>=perc].= perc;
    push!(post,aux)
    aux = round.(UInt8,255*(aux./maximum(aux)))
    push!(post2,aux)
end

## Plot
plot_cine(frames,floor(N_phases/1);Δt=TR)
plot_cine(post,floor(N_phases/1);Δt=TR)
plot_cine(post2,floor(N_phases/1);Δt=TR)


plot_cine(frames,N_phases) 

## Frame intensity profiling
tissue_profile = [frame[16,10]/1.6 for frame in post]
blood_profile  = [frame[16,16] for frame in post]
x = TR*(1:10)*1e3
blood = KomaMRIPlots.PlotlyJS.scatter(;x=x,y=blood_profile,name="Flow pixel intensity")
tissue = KomaMRIPlots.PlotlyJS.scatter(;x=x,y=tissue_profile,name="Wall tissue pixel intensity")


KomaMRIPlots.PlotlyJS.plot([tissue,blood])