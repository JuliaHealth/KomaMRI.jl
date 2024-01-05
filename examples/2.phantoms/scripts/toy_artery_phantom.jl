"""
This file generates a .phantom of a contracting and expanding ring.
A ventricle of the heart is recreated, with a diameter equal to FOV 
and an expansion and contraction of Δ meters every T seconds.
"""

using KomaMRI
using KomaMRICore.HDF5

FOV = 2e-2   # [m] Artery diameter
L = 4e-2     # [m] Artery length 
Δ = 0.2e-2   # [m] Displacement
T = [1]      # [s] Period
K = 16

#PARAMETERS
N = 21
Δxr = FOV/(N-1) #Aprox rec resolution, use Δx_pix and Δy_pix
Ns = 2 #number of spins per voxel
Δx = Δxr/sqrt(Ns) #spin separation
#POSITIONS
FOVx = FOVy = FOV
FOVz = L

x = -FOVx/2:Δx:FOVx/2-Δx 
y = -FOVy/2:Δx:FOVy/2-Δx 
z = -FOVz/2:Δx:FOVz/2-Δx 

xx = reshape(x, (length(x),1,1)) 
yy = reshape(y, (1,length(y),1)) 
zz = reshape(z, (1,1,length(z))) 

# Grid
x = 1*xx .+ 0*yy .+ 0*zz
y = 0*xx .+ 1*yy .+ 0*zz
z = 0*xx .+ 0*yy .+ 1*zz

#PHANTOM
⚪(R) =  (x.^2 .+ y.^2 .<= R^2) #Circle of radius R
v = FOV/4 #m/s 1/16 th of the FOV during acquisition
# Water spins
R = FOV/2
r = 6/11*FOV/2

tissue = ⚪(R) 
blood  = ⚪(r)

ρ = 0.9*Int.(tissue .| blood) #proton density

T1 = (1026*Int.(tissue .| blood)) * 1e-3   
T2 = (42*Int.(tissue .| blood))   * 1e-3 

# T1 = (1026*tissue + 1200*blood) * 1e-3   
# T2 = (42*tissue  + 100*blood)   * 1e-3 

x_values  = x[ρ .!= 0]
y_values  = y[ρ .!= 0]
z_values  = z[ρ .!= 0]

T1_values = T1[ρ .!= 0]
T2_values = T2[ρ .!= 0]
ρ_values  = ρ[ρ .!= 0]


# Blood flow
Δ = (maximum(z_values)-minimum(z_values))/K
Δ = ones(length(z_values),K-1)*Δ
Δ = cumsum(Δ,dims=2)

zt  = repeat(z_values,1,K-1)
ztt = zt + Δ
ztt[ztt .> maximum(z_values)] .-= (maximum(z_values)-minimum(z_values))
Δz = ztt - zt
Δz[x_values.^2 .+ y_values.^2 .> r^2,:] .= zeros(1,K-1)


reset_mag = BitMatrix(zeros(length(z_values),K))
for i in (1:size(Δz)[1])
    row = Δz[i,:]
    idx = findfirst(x -> x < 0, row)
    if isnothing(idx)
        reset_mag[i,end] = 1
    else
        if sum(abs.(row)) != 0
            reset_mag[i,idx] = 1
        end
    end
end

artery_vertical = Phantom(name="Flow Artery (Vertical)",
                          x=x_values,
                          y=y_values,
                          z=z_values,
                          T1=T1_values,
                          T2=T2_values,
                          ρ=ρ_values)

artery_vertical.mov = ArbitraryMotion(dur=[1.0],K=K,Δx=zeros(length(x_values),K-1),Δz=Δz,resetmag=reset_mag)


artery_horizontal = Phantom(name="Flow Artery (Horizontal)",
                            x=z_values,
                            y=y_values,
                            z=x_values,
                            T1=T1_values,
                            T2=T2_values,
                            ρ=ρ_values)

artery_horizontal.mov = ArbitraryMotion(dur=[1.0],K=K,Δx=Δz,resetmag=reset_mag)


write_phantom(artery_vertical, String(@__DIR__)*"/../flow_artery_vertical.phantom");
write_phantom(artery_horizontal, String(@__DIR__)*"/../flow_artery_horizontal.phantom");


#=
# Create HDF5 phantom file
fid = h5open(String(@__DIR__)*"/../toy_artery.phantom","w")

# Root attributes
attributes(fid)["Version"] = "1.0"
attributes(fid)["Name"] = "Toy Artery"
attributes(fid)["Ns"] = length(x_values)
attributes(fid)["Dims"] = 3
attributes(fid)["Dynamic"] = 1     # 0=False, 1=True

# Spin initial positions
pos = create_group(fid,"position")
x = create_group(pos,"x"); 
y = create_group(pos,"y"); 
z = create_group(pos,"z");

x["values"] = x_values
y["values"] = y_values
z["values"] = z_values


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
attributes(segments)["N"] = length(T)
attributes(segments)["K"] = K
segments["dur"] = T

motionz = create_group(motion,"motionz")
attributes(motionz)["type"] = "Explicit"
motionz["values"] = Δz

resetmag = create_group(motion,"resetmag")
attributes(resetmag)["type"] = "Explicit"
resetmag["values"] = Int.(reset_mag)

close(fid)
=#
