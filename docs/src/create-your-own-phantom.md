# Create Your Own Phantom

!!! warning
    This section is currently under construction, and some details on how to construct a Phantom may be missing.

This is an example of how to create a **Phantom** struct starting from a `.mat` file:

```julia
# Export necessary modules
using KomaMRI, MAT

# Create the function that creates a phantom
function phantom_example(; axis="axial", ss=4)

    # Get data from .mat file
    path_koma = dirname(dirname(pathof(KomaMRI)))
    path_phantom_mat = joinpath(path_koma, "KomaMRICore", "src", "datatypes", "phantom", "brain2D.mat")
    data = MAT.matread(path_phantom_mat)
    class = data[axis][1:ss:end,1:ss:end]

    # Define spin position vectors
    Δx = .5e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    # Define spin property vectors
    T2 = (class.==23)*329 .+ #CSF
        (class.==46)*83 .+ #GM
        (class.==70)*70 .+ #WM
        (class.==93)*70 .+ #FAT1
        (class.==116)*47 .+ #MUSCLE
        (class.==139)*329 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*70 .+ #FAT2
        (class.==232)*329 .+ #DURA
        (class.==255)*70 #MARROW
    T2s = (class.==23)*58 .+ #CSF
        (class.==46)*69 .+ #GM
        (class.==70)*61 .+ #WM
        (class.==93)*58 .+ #FAT1
        (class.==116)*30 .+ #MUSCLE
        (class.==139)*58 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*61 .+ #FAT2
        (class.==232)*58 .+ #DURA
        (class.==255)*61 .+#MARROW
        (class.==255)*70 #MARROW
    T1 = (class.==23)*2569 .+ #CSF
        (class.==46)*833 .+ #GM
        (class.==70)*500 .+ #WM
        (class.==93)*350 .+ #FAT1
        (class.==116)*900 .+ #MUSCLE
        (class.==139)*569 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*500 .+ #FAT2
        (class.==232)*2569 .+ #DURA
        (class.==255)*500 #MARROW
    ρ = (class.==23)*1 .+ #CSF
        (class.==46)*.86 .+ #GM
        (class.==70)*.77 .+ #WM
        (class.==93)*1 .+ #FAT1
        (class.==116)*1 .+ #MUSCLE
        (class.==139)*.7 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*.77 .+ #FAT2
        (class.==232)*1 .+ #DURA
        (class.==255)*.77 #MARROW
	Δw_fat = -220*2π
	Δw = (class.==93)*Δw_fat .+ #FAT1
		(class.==209)*Δw_fat    #FAT2
	T1 = T1*1e-3
	T2 = T2*1e-3
	T2s = T2s*1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "brain2D_"*axis,
		x = y[ρ.!=0],
		y = x[ρ.!=0],
		z = 0*x[ρ.!=0],
		ρ = ρ[ρ.!=0],
		T1 = T1[ρ.!=0],
		T2 = T2[ρ.!=0],
		T2s = T2s[ρ.!=0],
		Δw = Δw[ρ.!=0],
    )
    return obj
end

# Call the function to create a phantom
obj = phantom_example()

# Plot the proton density of the phantom
plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/create-your-own-phantom.html" style="width:100%; height:620px;"></object>
```