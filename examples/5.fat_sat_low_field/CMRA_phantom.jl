B0 = 0.55         # 0.55 T
fat_ppm = -3.4e-6 # 3.4ppm
Niso = 200        # 200 isochromats
Δx_voxel = 1.5e-3 # 1.5 mm
fat_freq = γ*B0*fat_ppm
dx = Array(range(-Δx_voxel/2, Δx_voxel/2, Niso))
#Phantom, values from: 
# Opportunities in Interventional and Diagnostic Imaging by Using High-Performance Low-Field-Strength MRI
# Adrienne E. Campbell-Washburn, Rajiv Ramasawmy, Matthew C. Restivo, Ipshita Bhattacharya, Burcu Basar, Daniel A. Herzka, Michael S. Hansen, Toby Rogers, W. Patricia Bandettini, Delaney R. McGuirt, Christine Mancini, David Grodzki, Rainer Schneider, Waqas Majeed, Himanshu Bhat, Hui Xue, Joel Moss, Ashkan A. Malayeri, Elizabeth C. Jones, Alan P. Koretsky, Peter Kellman, Marcus Y. Chen, Robert J. Lederman, and Robert S. Balaban
# Radiology 2019 293:2, 384-393
function cardiac_phantom(off; off_fat=fat_freq)
    myocard = Phantom{Float64}(x=dx, ρ=0.6*ones(Niso), T1=701e-3*ones(Niso),  
                               T2=58e-3*ones(Niso),    Δw=2π*off*ones(Niso))
    blood =   Phantom{Float64}(x=dx, ρ=0.7*ones(Niso), T1=1122e-3*ones(Niso), 
                               T2=263e-3*ones(Niso),   Δw=2π*off*ones(Niso))
    fat1 =    Phantom{Float64}(x=dx, ρ=1.0*ones(Niso), T1=183e-3*ones(Niso),
                               T2=93e-3*ones(Niso),    Δw=2π*(off_fat + off)*ones(Niso))
    fat2 =    Phantom{Float64}(x=dx, ρ=1.0*ones(Niso), T1=130e-3*ones(Niso),
                               T2=93e-3*ones(Niso),    Δw=2π*(off_fat + off)*ones(Niso))
    obj = myocard + blood + fat1 + fat2
    return obj
end