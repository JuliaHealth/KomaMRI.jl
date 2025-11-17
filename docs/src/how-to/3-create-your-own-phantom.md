# Create Your Own Phantom

In this section, we will create some custom **Phantom** structs.

In **KomaMRI**, the creation of a **Phantom** struct involves defining spin position arrays (x, y, z) and spin property arrays. 
The indices of these arrays are then associated with independent spins. See the [Phantom explanation section](../explanation/lit-1-phantom.jl) for more
information about how digital phantoms work in KomaMRI. 

## Basic case: 1-spin phantom

You can create a **Phantom** with one spin like so:
```julia
# Define arrays of positions (spin at zero position)
x = [0.0]
y = [0.0]
z = [0.0]

# Define arrays of properties (for CSF tissue)
ρ = [1.0]
T1 = [2.569]
T2 = [0.329]
T2s = [0.058]

# Define the phantom
spin = Phantom(name="spin", x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s)
```
```julia-repl
Phantom{Float64}
  name: String "spin"
  x: Array{Float64}((1,)) [0.0]
  y: Array{Float64}((1,)) [0.0]
  z: Array{Float64}((1,)) [0.0]
  ρ: Array{Float64}((1,)) [1.0]
  T1: Array{Float64}((1,)) [2.569]
  T2: Array{Float64}((1,)) [0.329]
  T2s: Array{Float64}((1,)) [0.058]
  Δw: Array{Float64}((1,)) [0.0]
  Dλ1: Array{Float64}((1,)) [0.0]
  Dλ2: Array{Float64}((1,)) [0.0]
  Dθ: Array{Float64}((1,)) [0.0]
  motion: NoMotion NoMotion()
```

You can add more properties to the **Phantom**, such as off-resonance, diffusion parameters, and even [motion](../explanation/2-motion.md) information.

## Loading phantom data from a .mat file
If you are familiar with the **MRI** world, you likely have a 2D or 3D array, where each element contains an ID number identifying a different class of tissue. In this setup, the array axes represent spatial positions, while the elements are used for tissue identification.

In this example, we will utilize a `.mat` file containing arrays with such arrangements. The file is readily available upon installing **KomaMRI**. Let's read the file and store the 2D data in an array called `class`:

```julia
# Get data from a .mat file
path_koma = dirname(dirname(pathof(KomaMRI)))
path_phantom_mat = joinpath(path_koma, "KomaMRIBase", "src", "datatypes", "phantom", "pelvis2D.mat")
data = KomaMRIFiles.matread(path_phantom_mat)
class = data["pelvis3D_slice"]
```

You can visualize the tissue map using the [`plot_image`](@ref) function:
```julia
plot_image(class)
```
```@raw html
<center><object type="text/html" data="../../assets/create-your-own-phantom-class-map.html" style="width:85%; height:470px;"></object></center>
```

Let's define the position arrays. You need to know the distance between the spins in the original array (in this case, it is 0.5mm), and then you can determine all the positions like this (the z-component is not calculated since this is a 2D example):
```julia
# Define spin position arrays
Δx = .5e-3                  # 0.5mm
M, N = size(class)          # Number of spins in x and y
FOVx = (M-1)*Δx             # Field of view in x
FOVy = (N-1)*Δx             # Field of view in y
x = -FOVx/2:Δx:FOVx/2       # x spin coordinates vector
y = -FOVy/2:Δx:FOVy/2       # y spin coordinates vector
x, y = x .+ y'*0, x*0 .+ y' # x and y grid points
```

Now, let's define the arrays for the properties. It's essential to have prior knowledge of the property values for different tissue classes. For example, for soft tissue, we use `ρ = 0.9`, `T1 = 1200 * 1e-3`, `T2 = 80 * 1e-3`, and `T2s = 80 * 1e-3`. Additionally, we create an array mask to identify the location of a tissue's ID. For soft tissue with ID = 153, the mask is `(class .== 153)`. Finally, to obtain a property, sum all the masks with values for all tissue classes. This process is illustrated below: 
```julia
# Define the proton density array
ρ = (class.==102)*.86 .+    # Fat
    (class.==153)*.9 .+     # SoftTissue
    (class.==204)*.4 .+     # SpongyBone
    (class.==255)*.2        # CorticalBone

# Define the T1 decay array
T1 = (class.==102)*366 .+   # Fat
    (class.==153)*1200 .+   # SoftTissue
    (class.==204)*381 .+    # SpongyBone
    (class.==255)*100       # CorticalBone

# Define the T2 decay array
T2 = (class.==102)*70 .+    # Fat
    (class.==153)*80 .+     # SoftTissue
    (class.==204)*52 .+     # SpongyBone
    (class.==255)*.3        # CorticalBone

# Define the T2s decay array
T2s = (class.==102)*70 .+   # Fat
    (class.==153)*80 .+     # SoftTissue
    (class.==204)*52 .+     # SpongyBone
    (class.==255)*.3        # CorticalBone

# Define off-resonance array
Δw_fat = -220 * 2π
Δw = (class.==102) * Δw_fat # FAT1

# Adjust with scaling factor
T1 = T1*1e-3
T2 = T2*1e-3
T2s = T2s*1e-3
```

Finally, we can invoke the [`Phantom`](@ref KomaMRIBase.Phantom) constructor. However, before doing so, we choose not to store spins where the proton density is zero to avoid unnecessary data storage. This is achieved by applying the mask `ρ.!=0` to the arrays. Additionally, please note that we set the z-position array filled with zeros.
```julia
# Define the phantom
obj = Phantom(
    name = "custom-pelvis",
	x = x[ρ.!=0],
	y = y[ρ.!=0],
	z = 0*x[ρ.!=0],
	ρ = ρ[ρ.!=0],
	T1 = T1[ρ.!=0],
	T2 = T2[ρ.!=0],
	T2s = T2s[ρ.!=0],
	Δw = Δw[ρ.!=0],
)
```

We can display the **Phantom** struct with the [`plot_phantom_map`](@ref) function. In this case we select the T1 decay to be displayed, but you can choose other property to be displayed:
```julia
plot_phantom_map(obj, :T1)
```
```@raw html
<object type="text/html" data="../../assets/create-your-own-phantom-pelvis-T1.html" style="width:85%; height:470px;"></object>
```

## Creating a custom flow cylinder phantom

This section shows how to create a slightly more complex phantom, which includes flow, from scratch. 

The phantom is modelled as a vertical cylindrical tube, with outer radius of 10 mm, an inner radius of 4.5 mm, and a length of 40 mm.

First, we generate the 3D grid that defines the spatial layaout of the spins:

```julia
R = 10e-3; r = 4.5e-3; L = 40e-3; Δx = 1e-3

#POSITIONS
x = -R:Δx:R
y = -R:Δx:R 
z = -L/2:Δx:L/2

xx = reshape(x, (length(x),1,1)) 
yy = reshape(y, (1,length(y),1)) 
zz = reshape(z, (1,1,length(z))) 

# Grid
x = 1*xx .+ 0*yy .+ 0*zz
y = 0*xx .+ 1*yy .+ 0*zz
z = 0*xx .+ 0*yy .+ 1*zz
```

Then, we separate the grid into two regions. The boolean array `bl` marks the inner part of the tube, where the "blood" spins should flow, 
while the `ts` array identifies the surrounding tissue, which must remain static. The `⚪(R)` function helps us with this task:

```julia
⚪(R) =  (x.^2 .+ y.^2 .<= R^2) # circle of radius R

ts = Bool.(⚪(R) - ⚪(r))
bl = Bool.(⚪(r))
```

Once each region is defined, we can go ahead and create the phantom. It is often convenient to generate the phantoms for each region separately. Let's start with the static outer `tissue`, which is very straightforward. We’ll assign it PD, T1, and T2 values of 1.0, 700 ms, and 42 ms, respectively:

```julia
PD = 1.0; T1 = 700e-3; T2 = 42e-3

tissue = Phantom(
    name="Tissue",
    x=x[ts],
    y=y[ts],
    z=z[ts],
    ρ=PD.*ones(length(x[ts])),
    T1=T1.*ones(length(x[ts])),
    T2=T2.*ones(length(x[ts]))
)
```

```julia-repl
Phantom{Float64}
  name: String "Tissue"
  x: Array{Float64}((10168,)) [0.0, -0.004, -0.003, -0.002, -0.001  …  0.0, 0.001, 0.002, 0.003, 0.004, 0.0]
  y: Array{Float64}((10168,)) [-0.01, -0.009, -0.009, -0.009, -0.009  …  0.009, 0.009, 0.009, 0.009, 0.01]
  z: Array{Float64}((10168,)) [-0.02, -0.02, -0.02, -0.02, -0.02  …  0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
  ρ: Array{Float64}((10168,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  …  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  T1: Array{Float64}((10168,)) [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7  …  0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7]
  T2: Array{Float64}((10168,)) [0.042, 0.042, 0.042, 0.042, 0.042  …  0.042, 0.042, 0.042, 0.042, 0.042]
  T2s: Array{Float64}((10168,)) [1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6  …  1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]
  Δw: Array{Float64}((10168,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dλ1: Array{Float64}((10168,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dλ2: Array{Float64}((10168,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dθ: Array{Float64}((10168,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  motion: NoMotion NoMotion()
```

Now, we repeat the process for the `blood` spins inside the tube. For this region, the PD, T1, and T2 values will be 0.9, 1200 ms, and 92 ms, respectively.

To define the flow movement for these inner spins, we need to create a [`FlowPath`](@ref) motion. This will include both the trajectories followed by each spin and the flags that indicate when a spin leaves the phantom's domain (at the upper end) and is reinjected at the lower end of the cylinder:

```julia
PD = 0.9; T1 = 1200e-3; T2 = 92e-3

Nt = 500    # Number of discrete time nodes
v  = 200e-3 # 200 mm/s 

dx = dy = zeros(length(z[bl]), Nt)
dz =z[bl] .+ cumsum(L/Nt .+ zeros(1,Nt), dims=2)

spin_reset = dz .> L/2
for i in 1:size(spin_reset, 1)
    idx = findfirst(x -> x == 1, spin_reset[i, :])
    if idx !== nothing
        spin_reset[i, :]  .= 0
        spin_reset[i, idx] = 1 # Activate the flag at the NEXT node after the jump
    end
end

dz[dz .> L/2] .-= L
dz .-= z[bl]

blood = Phantom(
    name="Blood",
    x=x[bl],
    y=y[bl],
    z=z[bl],
    ρ =PD.*ones(length(x[bl])),
    T1=T1.*ones(length(x[bl])),
    T2=T2.*ones(length(x[bl])),
    motion=FlowPath(dx, dy, dz, spin_reset, Periodic(L/v, 1.0-1e-6))
)
```

```julia-repl
Phantom{Float64}
  name: String "Blood"
  x: Array{Float64}((2829,)) [-0.002, -0.001, 0.0, 0.001, 0.002  …  -0.002, -0.001, 0.0, 0.001, 0.002]
  y: Array{Float64}((2829,)) [-0.004, -0.004, -0.004, -0.004  …  0.004, 0.004, 0.004, 0.004, 0.004]
  z: Array{Float64}((2829,)) [-0.02, -0.02, -0.02, -0.02, -0.02  …  0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
  ρ: Array{Float64}((2829,)) [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9  …  0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
  T1: Array{Float64}((2829,)) [1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2  …  1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
  T2: Array{Float64}((2829,)) [0.092, 0.092, 0.092, 0.092, 0.092  …  0.092, 0.092, 0.092, 0.092, 0.092]
  T2s: Array{Float64}((2829,)) [1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6  …  1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e6]
  Δw: Array{Float64}((2829,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dλ1: Array{Float64}((2829,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dλ2: Array{Float64}((2829,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  Dθ: Array{Float64}((2829,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  motion: Motion{Float64}
```

Last, we obtain the final `obj` as the sum of both phantoms:

```julia
obj = tissue + blood

plot_phantom_map(obj, :T1; time_samples=11)
```
```@raw html
<object type="text/html" data="../../assets/create-your-own-phantom-flow-T1.html" style="width:85%; height:470px;"></object>
```

## Importing/exporting phantoms (.phantom file format)
Another option to create your own phantom is to import it directly from a [`.phantom`](../explanation/3-phantom-format.md) file.
To demonstrate this, we will first create the file from the previous phantom, using the [`write_phantom`](@ref) function, and then we will read it using the [`read_phantom`](@ref) function:

```julia-repl
julia> write_phantom(obj, "flow_cylinder.phantom")
# Now, you could even close Julia, turn off your computer...
julia> obj2 = read_phantom("flow_cylinder.phantom")
```