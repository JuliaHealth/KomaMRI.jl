# Create Your Own Phantom

In this section, we will create a custom **Phantom** struct. While the example is presented in 2D, the concepts discussed here can be readily extended to 3D phantoms.

In **KomaMRI**, the creation of a **Phantom** struct involves defining spin position arrays (x, y, z) and spin property arrays. The indices of these arrays are then associated with independent spins.

For instance, you can create a **Phantom** with one spin like so:
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
  ux: #122 (function of type KomaMRICore.var"#122#136")
  uy: #123 (function of type KomaMRICore.var"#123#137")
  uz: #124 (function of type KomaMRICore.var"#124#138")
```

You can add more properties to the **Phantom**, such as off-resonance, diffusion parameters, and even functions of motion. However, we won't be utilizing them (except for the off-resonance parameter) to maintain simplicity.

If you are familiar with the **MRI** world, you likely have a 2D or 3D array, where each element contains an ID number identifying a different class of tissue. In this setup, the array axes represent spatial positions, while the elements are used for tissue identification.

In this example, we will utilize a `.mat` file containing arrays with such arrangements. The file is readily available upon installing **KomaMRI**. Let's read the file and store the 2D data in an array called `class`:"
```julia
# Import necessary modules
using KomaMRI, MAT

# Get data from a .mat file
path_koma = dirname(dirname(pathof(KomaMRI)))
path_phantom_mat = joinpath(path_koma, "KomaMRICore", "src", "datatypes","phantom", "brain2D.mat")
data = MAT.matread(path_phantom_mat)
class = data["sagittal"]
```

You can visualize the tissue map using the [`plot_image`](@ref) function:
```julia
plot_image(class)
```
```@raw html
<center><object type="text/html" data="../assets/create-your-own-phantom-class-map.html" style="width:100%; height:620px;"></object></center>
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

Now, let's define the arrays for the properties. It's essential to have prior knowledge of the property values for different tissue classes. For example, for muscle tissue, we use `ρ = 1`, `T1 = 900 * 1e-3`, `T2 = 47 * 1e-3`, and `T2s = 30 * 1e-3`. Additionally, create an array mask to identify the location of a tissue's ID. For muscle with ID = 116, the mask is `(class == 116)`. Finally, to obtain a property, sum all the masks with values for all tissue classes. This process is illustrated below: 
```julia
# Define the proton density array
ρ = (class.==23)*1 .+       # CSF
    (class.==46)*.86 .+     # GM
    (class.==70)*.77 .+     # WM
    (class.==93)*1 .+       # FAT1
    (class.==116)*1 .+      # MUSCLE
    (class.==139)*.7 .+     # SKIN/MUSCLE
    (class.==162)*0 .+      # SKULL
    (class.==185)*0 .+      # VESSELS
    (class.==209)*.77 .+    # FAT2
    (class.==232)*1 .+      # DURA
    (class.==255)*.77       # MARROW

# Define the T1 decay array
T1 = (class.==23)*2569 .+   # CSF
    (class.==46)*833 .+     # GM
    (class.==70)*500 .+     # WM
    (class.==93)*350 .+     # FAT1
    (class.==116)*900 .+    # MUSCLE
    (class.==139)*569 .+    # SKIN/MUSCLE
    (class.==162)*0 .+      # SKULL
    (class.==185)*0 .+      # VESSELS
    (class.==209)*500 .+    # FAT2
    (class.==232)*2569 .+   # DURA
    (class.==255)*500       # MARROW

# Define the T2 decay array
T2 = (class.==23)*329 .+    # CSF
    (class.==46)*83 .+      # GM
    (class.==70)*70 .+      # WM
    (class.==93)*70 .+      # FAT1
    (class.==116)*47 .+     # MUSCLE
    (class.==139)*329 .+    # SKIN/MUSCLE
    (class.==162)*0 .+      # SKULL
    (class.==185)*0 .+      # VESSELS
    (class.==209)*70 .+     # FAT2
    (class.==232)*329 .+    # DURA
    (class.==255)*70        # MARROW

# Define the T2s decay array
T2s = (class.==23)*58 .+    # CSF
    (class.==46)*69 .+      # GM
    (class.==70)*61 .+      # WM
    (class.==93)*58 .+      # FAT1
    (class.==116)*30 .+     # MUSCLE
    (class.==139)*58 .+     # SKIN/MUSCLE
    (class.==162)*0 .+      # SKULL
    (class.==185)*0 .+      # VESSELS
    (class.==209)*61 .+     # FAT2
    (class.==232)*58 .+     # DURA
    (class.==255)*61        # MARROW

# Define off-resonance array
Δw_fat = -220*2π
Δw = (class.==93)*Δw_fat .+ # FAT1
    (class.==209)*Δw_fat    # FAT2

# Adjust with scaling factor
T1 = T1*1e-3
T2 = T2*1e-3
T2s = T2s*1e-3
```

Finally, we can invoke the [`Phantom`](@Ref) constructor. However, before doing so, we choose not to store spins where the proton density is zero to avoid unnecessary data storage. This is achieved by applying the mask `ρ.!=0` to the arrays. Additionally, please note that we set the x-position array filled with zeros, and we interchange the y and z coordinates.
```julia
# Define the phantom
obj = Phantom{Float64}(
    name = "custom-brain",
	x = 0*x[ρ.!=0],
	y = x[ρ.!=0],
	z = y[ρ.!=0],
	ρ = ρ[ρ.!=0],
	T1 = T1[ρ.!=0],
	T2 = T2[ρ.!=0],
	T2s = T2s[ρ.!=0],
	Δw = Δw[ρ.!=0],
)
```

We can display the **Phantom** struct with the [`plot_phantom_map`](@ref) function. In this case we select the proton density to be displayed, but you can choose other property to be displayed:
```julia
plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/create-your-own-phantom-plot-rho.html" style="width:100%; height:620px;"></object>
```
