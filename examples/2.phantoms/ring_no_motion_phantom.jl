"""
This file generates a .phantom of an static ring.
"""

using HDF5

#PARAMETERS
FOV = 10e-2 #m Diameter ventricule
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

# Create HDF5 phantom file
fid = h5open(String(@__DIR__)*"/ring_no_motion.phantom","w")

# Root attributes
attributes(fid)["Version"] = "1.0"
attributes(fid)["Name"] = "Ring No Motion"
attributes(fid)["Ns"] = length(x_values)
attributes(fid)["Dims"] = 2
attributes(fid)["Dynamic"] = 0     # 0=False, 1=True

# Spin initial positions
pos = create_group(fid,"position")
x = create_group(pos,"x"); 
y = create_group(pos,"y"); 

attributes(x)["type"] = "Explicit"
attributes(y)["type"] = "Explicit"

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

close(fid)
