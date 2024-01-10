"""
    obj = Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, ux, uy, uz)

The `Phantom` struct has field names that are vectors associated with spin properties. It
allows accessing specific spin property values by indexing these vectors. This struct serves
as input for the simulation.

# Arguments
- `name`: (`::String`) phantom name
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `ρ`: (`::AbstractVector{T<:Real}`) spin proton density vector
- `T1`: (`::AbstractVector{T<:Real}`, `[s]`) spin T1 parameter vector
- `T2`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2 parameter vector
- `T2s`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2s parameter vector
- `Δw`: (`::AbstractVector{T<:Real}`, `[rad/s]`) spin off-resonance parameter vector
- `Dλ1`: (`::AbstractVector{T<:Real}`) spin Dλ1 (diffusion) parameter vector
- `Dλ2`: (`::AbstractVector{T<:Real}`) spin Dλ2 (diffusion) parameter vector
- `Dθ`: (`::AbstractVector{T<:Real}`) spin Dθ (diffusion) parameter vector
- `ux`: (`::Function`) displacement field in the x-axis
- `uy`: (`::Function`) displacement field in the y-axis
- `uz`: (`::Function`) displacement field in the z-axis

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = Phantom(x=[0.0])

julia> obj.ρ
```
"""
 @with_kw mutable struct Phantom{T<:Real}
    name::String = "spins"
	x::AbstractVector{T}
	y::AbstractVector{T} = zeros(size(x))
	z::AbstractVector{T} = zeros(size(x))
	ρ::AbstractVector{T} = ones(size(x))
	T1::AbstractVector{T} = ones(size(x)) * 1_000_000
	T2::AbstractVector{T} = ones(size(x)) * 1_000_000
	T2s::AbstractVector{T} = ones(size(x)) * 1_000_000
	#Off-resonance related
	Δw::AbstractVector{T} = zeros(size(x))
	#χ::Vector{SusceptibilityModel}
	#Diffusion
	Dλ1::AbstractVector{T} = zeros(size(x))
	Dλ2::AbstractVector{T} = zeros(size(x))
	Dθ::AbstractVector{T} =  zeros(size(x))
	#Diff::Vector{DiffusionModel}  #Diffusion map
	#Motion
	ux::Function = (x,y,z,t)->0
	uy::Function = (x,y,z,t)->0
	uz::Function = (x,y,z,t)->0
end

"""Size and length of a phantom"""
size(x::Phantom) = size(x.ρ)
Base.length(x::Phantom) = length(x.ρ)
# To enable to iterate and broadcast over the Phantom
Base.iterate(x::Phantom) = (x[1], 2)
Base.iterate(x::Phantom, i::Integer) = (i <= length(x)) ? (x[i], i+1) : nothing
Base.lastindex(x::Phantom) = length(x)
Base.getindex(x::Phantom, i::Integer) = x[i:i]

"""Compare two phantoms"""
Base.isapprox(obj1::Phantom, obj2::Phantom)  = begin
    obj1.x     ≈ obj2.x    &&
    obj1.y     ≈ obj2.y    &&
    obj1.z     ≈ obj2.z    &&
    obj1.ρ     ≈ obj2.ρ    &&
    obj1.T1    ≈ obj2.T1   &&
    obj1.T2    ≈ obj2.T2   &&
    obj1.T2s   ≈ obj2.T2s  &&
    obj1.Δw    ≈ obj2.Δw   &&
    obj1.Dλ1   ≈ obj2.Dλ1  &&
    obj1.Dλ2   ≈ obj2.Dλ2  &&
    obj1.Dθ    ≈ obj2.Dθ
end

"""
Separate object spins in a sub-group
"""
Base.getindex(obj::Phantom, p::AbstractRange) = begin
	Phantom(name=obj.name,
			x=obj.x[p],
			y=obj.y[p],
			z=obj.z[p],
			ρ=obj.ρ[p],
			T1=obj.T1[p],
			T2=obj.T2[p],
			T2s=obj.T2s[p],
			Δw=obj.Δw[p],
			#Diff=obj.Diff[p], #TODO!
			Dλ1=obj.Dλ1[p],
			Dλ2=obj.Dλ2[p],
			Dθ=obj.Dθ[p],
			#Χ=obj.Χ[p], #TODO!
			ux=obj.ux,
			uy=obj.uy,
			uz=obj.uz
			)
end

"""Separate object spins in a sub-group (lightweigth)."""
Base.view(obj::Phantom, p::AbstractRange) = begin
	@views Phantom(name=obj.name,
			x=obj.x[p],
			y=obj.y[p],
			z=obj.z[p],
			ρ=obj.ρ[p],
			T1=obj.T1[p],
			T2=obj.T2[p],
			T2s=obj.T2s[p],
			Δw=obj.Δw[p],
			#Diff=obj.Diff[p], #TODO!
			Dλ1=obj.Dλ1[p],
			Dλ2=obj.Dλ2[p],
			Dθ=obj.Dθ[p],
			#Χ=obj.Χ[p], #TODO!
			ux=obj.ux,
			uy=obj.uy,
			uz=obj.uz
			)
end

"""Addition of phantoms"""
+(s1::Phantom,s2::Phantom) = begin
	Phantom(name=s1.name*"+"*s2.name,
		x=[s1.x;s2.x],
		y=[s1.y;s2.y],
		z=[s1.z;s2.z],
		ρ=[s1.ρ;s2.ρ],
		T1=[s1.T1;s2.T1],
		T2=[s1.T2;s2.T2],
		T2s=[s1.T2s;s2.T2s],
		Δw=[s1.Δw;s2.Δw],
		#Diff=obj.Diff[p], #TODO!
		Dλ1=[s1.Dλ1;s2.Dλ1],
		Dλ2=[s1.Dλ2;s2.Dλ2],
		Dθ=[s1.Dθ;s2.Dθ],
		#Χ=obj.Χ[p], #TODO!
		ux=s1.ux,
		uy=s1.uy,
		uz=s1.uz
	)
end

"""Scalar multiplication of a phantom"""
*(α::Real,obj::Phantom) = begin
	Phantom(name=obj.name,
		x=obj.x,
		y=obj.y,
		z=obj.z,
		ρ=α*obj.ρ, #Only affects the proton density
		T1=obj.T1,
		T2=obj.T2,
		T2s=obj.T2s,
		Δw=obj.Δw,
		#Diff=obj.Diff[p], #TODO!
		Dλ1=obj.Dλ1,
		Dλ2=obj.Dλ2,
		Dθ=obj.Dθ,
		#Χ=obj.Χ[p], #TODO!
		ux=obj.ux,
		uy=obj.uy,
		uz=obj.uz
	)
end

"""
    obj = brain_phantom2D(; axis="axial", ss=4)

Creates a two-dimensional brain Phantom struct.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb

# Keywords
- `axis`: (`::String`, `="axial"`, opts=[`"axial"`, `"coronal"`, `"sagittal"`]) orientation
    of the phantom
- `ss`: (`::Integer`, `=4`) subsampling parameter in all axis

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom2D(; axis="sagittal", ss=1)

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom2D(; axis="axial", ss=4)

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path*"/phantom/brain2D.mat")
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

"""
    obj = brain_phantom3D(; ss=4, start_end=[160, 200])

Creates a three-dimentional brain Phantom struct.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb

# Keywords
- `ss`: (`::Integer`, `=4`) subsampling parameter in all axes
- `start_end`: (`::Integer`, `=[160, 200]`) start and end indices along the z-axis.
    Acceptable values range from 1 to 362

# Returns
- `obj`: (`::Phantom`) 3D Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom3D(; ss=5)

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom3D(; ss=4, start_end=[160, 200])

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path*"/phantom/brain3D.mat")
    class = data["data"][1:ss:end,1:ss:end,start_end[1]:ss:start_end[2]]

    # Define spin position vectors
    Δx = .5e-3*ss
    M, N, Z = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
	FOVz = (Z-1)*Δx #[m]
    xx = reshape(-FOVx/2:Δx:FOVx/2,M,1,1) #spin coordinates
    yy = reshape(-FOVy/2:Δx:FOVy/2,1,N,1) #spin coordinates
	zz = reshape(-FOVz/2:Δx:FOVz/2,1,1,Z) #spin coordinates
    x = 1*xx .+ 0*yy .+ 0*zz
	y = 0*xx .+ 1*yy .+ 0*zz
	z = 0*xx .+ 0*yy .+ 1*zz

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
        name = "brain3D",
		x = y[ρ.!=0],
		y = x[ρ.!=0],
		z = z[ρ.!=0],
		ρ = ρ[ρ.!=0],
		T1 = T1[ρ.!=0],
		T2 = T2[ρ.!=0],
		T2s = T2s[ρ.!=0],
		Δw = Δw[ρ.!=0],
    )
	return obj
end

"""
    obj = pelvis_phantom2D(; ss=4)

Creates a two-dimensional pelvis Phantom struct.

# Keywords
- `ss`: (`::Integer`, `=4`) subsampling parameter

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = pelvis_phantom2D(; ss=1)

julia> pelvis_phantom2D(obj, :ρ)
```
"""
function pelvis_phantom2D(; ss=4)

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path*"/phantom/pelvis2D.mat")
    class = data["pelvis3D_slice"][1:ss:end,1:ss:end]

    # Define spin position vectors
    Δx = .5e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx             # [m]
    FOVy = (N-1)*Δx             # [m]
    x = -FOVx/2:Δx:FOVx/2       # spin coordinates
    y = -FOVy/2:Δx:FOVy/2       # spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' # grid points

    # Define spin property vectors
    ρ = (class.==51)*.001 .+    # Air
        (class.==102)*.86 .+    # Fat
        (class.==153)*.9 .+     # SoftTissue
        (class.==204)*.4 .+     # SpongyBone
        (class.==255)*.2        # CorticalBone
    T1 = (class.==51)*.001 .+   # Air
        (class.==102)*366 .+    # Fat
        (class.==153)*1200 .+   # SoftTissue
        (class.==204)*381 .+    # SpongyBone
        (class.==255)*100       # CorticalBone
    T2 = (class.==51)*.001 .+   # Air
        (class.==102)*70 .+     # Fat
        (class.==153)*80 .+     # SoftTissue
        (class.==204)*52 .+     # SpongyBone
        (class.==255)*.3        # CorticalBone
    T2s = (class.==51)*.001 .+  # Air
        (class.==102)*70 .+     # Fat
        (class.==153)*80 .+     # SoftTissue
        (class.==204)*52 .+     # SpongyBone
        (class.==255)*.3        # CorticalBone
	Δw_fat = -220 * 2π
	Δw = (class.==102) * Δw_fat # FAT1
	T1 = T1*1e-3
	T2 = T2*1e-3
	T2s = T2s*1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "pelvis2D",
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
