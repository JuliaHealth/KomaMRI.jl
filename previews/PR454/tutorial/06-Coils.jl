using KomaMRI, MAT # hide
sys = Scanner() # hide

ph_file = joinpath(dirname(pathof(KomaMRI)), "../examples/2.phantoms/sphere_fields.mat")
sphere = matread(ph_file)
Δx = 1e-3                   # 5mm
N = size(sphere["b1m"])  # Number of spins
FOV = (N .- 1) .* Δx    # Field of view
xr = -FOV[1]/2:Δx:FOV[1]/2  # x spin coordinates vector
ρ = abs.(sphere["b1m"][:]) .> 0
x = [x for (x, y, z) in Iterators.product(xr, xr, xr)][ρ .!= 0]
y = [y for (x, y, z) in Iterators.product(xr, xr, xr)][ρ .!= 0]
z = [z for (x, y, z) in Iterators.product(xr, xr, xr)][ρ .!= 0]
coil_sens = 1.0 * sphere["b1m"][:][ρ .!= 0] ./ maximum(abs.(sphere["b1m"][:][ρ .!= 0]))
ρ = 1.0 * ρ[ρ .!= 0]
obj = Phantom(; x, y, z, ρ, coil_sens)
p1 = plot_phantom_map(obj, :coil_sens ; height=400, width=400, darkmode=true)

display(p1)

using Interpolations
coil_sens = sphere["b1m"] ./ maximum(abs.(sphere["b1m"][:]))
obj = brain_phantom2D()
obj.coil_sens .= LinearInterpolation((xr,xr,xr), coil_sens, extrapolation_bc=0).(obj.x, obj.y, obj.z)
p2 = plot_phantom_map(obj, :coil_sens ; height=400, width=400, darkmode=true)
display(p2)

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seq_file)

raw = simulate(obj, seq, sys)
acq = AcquisitionData(raw) # hide
acq.traj[1].circular = false # hide
Nx, Ny = raw.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image = reconstruction(acq, reconParams) # hide
slice_abs = abs.(image[:, :, 1]) # hide
p3 = plot_image(slice_abs; height=400) # hide
display(p3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
