"""
    obj = Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, ux, uy, uz)

The Phantom struct.

# Arguments
- `name`: (`::String`) name of the phantom
- `x`: (`::AbstractVector{T}`, `[m]`) vector of x-positions of the spins
- `y`: (`::AbstractVector{T}`, `[m]`) vector of y-positions of the spins
- `z`: (`::AbstractVector{T}`, `[m]`) vector of z-positions of the spins
- `ρ`: (`::AbstractVector{T}`) vector of proton density of the spins
- `T1`: (`::AbstractVector{T}`, `[s]`) vector of T1 parameters of the spins
- `T2`: (`::AbstractVector{T}`, `[s]`) vector of T2 parameters of the spins
- `T2s`: (`::AbstractVector{T}`, `[s]`) vector of T2s parameters of the spins
- `Δw`: (`::AbstractVector{T}`, `[rad/s]`) vector of off-resonance parameters of the spins
- `Dλ1`: (`::AbstractVector{T}`) vector of Dλ1 (diffusion) parameters of the spins
- `Dλ2`: (`::AbstractVector{T}`) vector of Dλ2 (diffusion) parameters of the spins
- `Dθ`: (`::AbstractVector{T}`) vector of Dθ (diffusion) parameters of the spins
- `ux`: (`::Function`) displacement field in the x-axis
- `uy`: (`::Function`) displacement field in the y-axis
- `uz`: (`::Function`) displacement field in the z-axis

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = Phantom()

julia> obj.ρ
```
"""
 @with_kw mutable struct Phantom{T<:Real}
    name::String = "spins"
	x::AbstractVector{T}
	y::AbstractVector{T} =   zeros(size(x))
	z::AbstractVector{T} =   zeros(size(x))
	ρ::AbstractVector{T} =   ones(size(x))
	T1::AbstractVector{T} =  ones(size(x)) * 1_000_000
	T2::AbstractVector{T} =  ones(size(x)) * 1_000_000
	T2s::AbstractVector{T} = ones(size(x)) * 1_000_000
	#Off-resonance related
	Δw::AbstractVector{T} =  zeros(size(x))
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

"""Separate object spins in a sub-group."""
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
    phantom = brain_phantom2D(;axis="axial", ss=4)

Creates a two-dimentional brain phantom struct.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb

# Keywords
- `axis`: (`::String`, `="axial"`, opts=[`"axial"`]) orientation of the phantom
- `ss`: (`::Real`, `=4`) subsampling parameter in all axis

# Returns
- `phantom`: (`::Phantom`) 2D Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom2D()

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom2D(;axis="axial", ss=4)
    path = @__DIR__
    data = MAT.matread(path*"/phantom/brain2D.mat")

    class = data[axis][1:ss:end,1:ss:end]
    Δx = .5e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

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
    phantom = Phantom{Float64}(name="brain2D_"*axis,
				x=y[ρ.!=0],
				y=x[ρ.!=0],
				z=0*x[ρ.!=0],
				ρ=ρ[ρ.!=0],
				T1=T1[ρ.!=0],
				T2=T2[ρ.!=0],
				T2s=T2s[ρ.!=0],
				Δw=Δw[ρ.!=0])
	phantom
end

"""
    phantom = brain_phantom3D(;ss=4)

Creates a three-dimentional brain phantom struct.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb

# Keywords
- `ss`: (`::Real`, `=4`) subsampling parameter in all axis

# Returns
- `phantom`: (`::Phantom`) 3D Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom3D()

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom3D(;ss=4,start_end=[160, 200])
    path = @__DIR__
    data = MAT.matread(path*"/phantom/brain3D.mat")

    class = data["data"][1:ss:end,1:ss:end,start_end[1]:ss:start_end[2]]
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
    phantom = Phantom{Float64}(name ="brain3D",
				x = y[ρ.!=0],
				y = x[ρ.!=0],
				z = z[ρ.!=0],
				ρ = ρ[ρ.!=0],
				T1 = T1[ρ.!=0],
				T2 = T2[ρ.!=0],
				T2s = T2s[ρ.!=0],
				Δw = Δw[ρ.!=0])
	phantom
end
