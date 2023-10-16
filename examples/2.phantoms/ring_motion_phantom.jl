"""
This file generates a .phantom of a contracting and expanding ring.
A ventricle of the heart is recreated, with a diameter equal to FOV 
and an expansion and contraction of Δ meters every T seconds.
"""

using KomaMRICore.HDF5

FOV = 10e-2 # [m] Diameter ventricule
Δ = 1e-2    # [m] Displacement
T = 1       # [s] Period

#PARAMETERS
N = 21
Δxr = FOV/(N-1) #Aprox rec resolution, use Δx_pix and Δy_pix
Ns = 100 #number of spins per voxel
Δx = Δxr/sqrt(Ns) #spin separation
#POSITIONS
x = y = -FOV/2:Δx:FOV/2-Δx #spin coordinates
x, y = x .+ y'*0, x*0 .+ y' #grid points
#PHANTOM
⚪(R) =  (x.^2 .+ y.^2 .<= R^2)*1. #Circle of radius R
v = FOV/4 #m/s 1/16 th of the FOV during acquisition
# Water spins
R = 9/10*FOV/2
r = 6/11*FOV/2
ring = ⚪(R) .- ⚪(r)

ρ = 0.9*ring #proton density
T1 = (1026*ring)*1e-3 #Myocardial T1
T2 = (42*ring)*1e-3 #T2 map [s]

x_values  = x[ρ .!= 0]
y_values  = y[ρ .!= 0]
T1_values = T1[ρ .!= 0]
T2_values = T2[ρ .!= 0]
ρ_values  = ρ[ρ .!= 0]

α = atan.((abs.(y_values))./(abs.(x_values))) 
α[(x_values.<=0).&(y_values.>=0)] =  π  .- α[(x_values.<=0).&(y_values.>=0)] 
α[(x_values.<=0).&(y_values.<0)] =  π  .+ α[(x_values.<=0).&(y_values.<0)]
α[(x_values.>0).&(y_values.<=0)] = 2*π .- α[(x_values.>0).&(y_values.<=0)]

# Create HDF5 phantom file
fid = h5open(String(@__DIR__)*"/ring_motion.phantom","w")

# Root attributes
attributes(fid)["Version"] = "1.0"
attributes(fid)["Name"] = "Ring Motion"
attributes(fid)["Ns"] = length(x_values)
attributes(fid)["Dims"] = 2
attributes(fid)["Dynamic"] = 1     # 0=False, 1=True

# Spin initial positions
pos = create_group(fid,"position")
x = create_group(pos,"x"); 
y = create_group(pos,"y"); 

x["values"] = x_values
y["values"] = y_values


# Contrast (Rho, T1, T2, Deltaw)
contrast = create_group(fid,"contrast")

rho = create_group(contrast,"ρ")
attributes(rho)["type"] = "Explicit"
rho["values"] = ρ_values 

T1 = create_group(contrast,"T1")
attributes(T1)["type"] = "Explicit"
T1["values"] = T1_values 

T2 = create_group(contrast,"T2")
attributes(T2)["type"] = "Explicit"
T2["values"] = T2_values 

Deltaw = create_group(contrast,"Δw")
attributes(Deltaw)["type"] = "Explicit"
Deltaw["values"] = zeros(length(x_values))


# Motion
motion = create_group(fid,"motion")
attributes(motion)["model"] = "Arbitrary"

segments = create_group(motion, "segments")
attributes(segments)["N"] = 1
attributes(segments)["K"] = 4
segments["dur"] = [T]

motionx = create_group(motion,"motionx")
attributes(motionx)["type"] = "Explicit"
motionx["values"] = hcat(Δ*cos.(α), zeros(length(x_values)), -Δ*cos.(α))
# motionx["values"] = hcat(Δ*ones(length(x_values)), zeros(length(x_values)), -Δ*ones(length(x_values)))

motiony = create_group(motion,"motiony")
attributes(motiony)["type"] = "Explicit"
motiony["values"] = hcat(Δ*sin.(α), zeros(length(x_values)), -Δ*sin.(α))
# motiony["values"] = hcat(Δ*ones(length(x_values)), zeros(length(x_values)), -Δ*ones(length(x_values)))

close(fid)
