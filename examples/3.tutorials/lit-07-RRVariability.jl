# # RR Variability in Cardiac Cine

using KomaMRI # hide

using CUDA

include("../1.sequences/koma_scripts/bSFFP_cine.jl")
include("utils/plot_cine.jl")

sys = Scanner() # hide

# In this tutorial, we will explore the effects of RR interval variability during a cardiac cine MRI acquisition, 
# and how these effects can be mitigated using proper prospective gating.

# We will begin by simulating a cardiac cine on a myocardial phantom with a constant RR interval.
# Let's call the `heart_phantom` function to create a ring-shaped phantom filled with blood, which resembles the left ventricle:
obj = heart_phantom()

# By default, this phantom exhibits periodic contraction and rotation, with a period of 1 second:
p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21) # hide
#md savefig(p1, "../assets/tut-6-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/tut-6-phantom.html" style="width:85%; height:470px;"></object></center>
#md # ```

# Now, we will create a bSSFP cine sequence with the following parameters:
##
RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 40          # image size = N x N
N_phases     = 40          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 20e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 40, adc_duration = adc_duration,
)

## Simulation # hide
raw1 = simulate(obj, seq, sys) # hide

## Reconstruction
frames = []
@info "Running reconstruction"
@time begin
    recParams = Dict{Symbol,Any}(:reco=>"direct")
    Nx = Ny = N_matrix
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = false

    acqData = AcquisitionData(raw1)

    _, ktraj = get_kspace(seq)

    for i in 1:N_phases
        acqAux = copy(acqData)
        range = reduce(vcat,[j*(N_matrix*N_phases).+((i-1)*N_matrix.+(1:N_matrix)) for j in 0:N_matrix-1])

        # Kdata
        acqAux.kdata[1] = reshape(acqAux.kdata[1][range],(N_matrix^2,1))

        # Traj
        acqAux.traj[1].circular = false

        acqAux.traj[1].nodes = transpose(ktraj[:, 1:2])[:,range]
        acqAux.traj[1].nodes = acqAux.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqAux.traj[1].nodes[:]))

        acqAux.traj[1].numProfiles = N_matrix
        acqAux.traj[1].times = acqAux.traj[1].times[range]

        # subsampleIndices
        acqAux.subsampleIndices[1] = acqAux.subsampleIndices[1][1:N_matrix^2]

        # Reconstruction
        aux = @timed reconstruction(acqAux, recParams)
        image  = reshape(aux.value.data,Nx,Ny,:)
        image_aux = abs.(image[:,:,1])

        push!(frames,image_aux)
    end
end

fps = 25
p2 = plot_cine(frames, fps; Δt=TR, filename="../assets/tut-7-frames.gif") # hide
#jl display(p2)

#md # ```@raw html
#md # <center><img src="../assets/tut-7-frames.gif" style="width:85%"></object></center>
#md # ```
##

# RRs = [1.0, 0.8, 1.0, 1.0] # hide
# obj.motion.motions[1].time.periods = RRs # hide
# obj.motion.motions[2].time.periods = RRs # hide

