"""
    phantom = read_phantom_jemris(filename)

Returns the Phantom struct from JEMRIS phantoms a file `.h5`.

# Arguments
- `filename`: (`::String`) the absolute or relative path of the phantom file `.h5`

# Returns
- `phantom`: (`::Phantom`) the phantom struct

# Examples
```julia-repl
julia> obj = read_phantom_jemris("examples/2.phantoms/brain.h5")
Phantom
  name: String "brain.h5"
  x: Array{Float64}((25841,)) [-0.0085, -0.0075  …  0.0035]
  y: Array{Float64}((25841,)) [-0.0985, -0.0985  …  0.1055]
  z: Array{Float64}((25841,)) [0.0, 0.0  …  0.0]
  ρ: Array{Float64}((25841,)) [1.0, 1.0  …  1.0]
  T1: Array{Float64}((25841,)) [2.569, 2.569  …  2.569]
  T2: Array{Float64}((25841,)) [0.329, 0.329  …  0.329]
  T2s: Array{Float64}((25841,)) [Inf, Inf  …  Inf]
  Δw: Array{Float64}((25841,)) [0.0, 0.0  …  0.0]
  Dλ1: Array{Float64}((25841,)) [0.0, 0.0  …  0.0]
  Dλ2: Array{Float64}((25841,)) [0.0, 0.0  …  0.0]
  Dθ: Array{Float64}((25841,)) [0.0, 0.0  …  0.0]
  ux: #161 (function of type KomaMRI.var"#161#162"{Int64})
  uy: #387 (function of type KomaMRI.var"#387#395")
  uz: #388 (function of type KomaMRI.var"#388#396")

julia> plot_phantom_map(obj, :ρ)
```
"""
function read_phantom_jemris(filename)
	# A(:,:,:,1)=Sample.M0;
	# I=find(Sample.T1); R1 =zeros(size(Sample.T1));  R1(I) =1./Sample.T1(I);
	# I=find(Sample.T2); R2 =zeros(size(Sample.T2));  R2(I) =1./Sample.T2(I);
	# I=find(Sample.T2S);R2S=zeros(size(Sample.T2S)); R2S(I)=1./Sample.T2S(I);
	# A(:,:,:,2)=R1;	#1/T1
	# A(:,:,:,3)=R2;	#1/T2
	# A(:,:,:,4)=R2S;	#1/T2s
	# A(:,:,:,5)=Sample.DB;
	fid = HDF5.h5open(filename)
	data = read(fid["sample/data"])
	Δx = read(fid["sample/resolution"]) * 1e-3 #[m]
	offset = read(fid["sample/offset"]) * 1e-3 #[m]
	mask = data[1,:,:,:] .!= 0
	#Maps
	ρ =   data[1,:,:,:]
	T1 =  1e-3 ./ data[2,:,:,:]
	T2 =  1e-3 ./ data[3,:,:,:]
	T2s = 1e-3 ./ data[4,:,:,:]
	Δw = data[5,:,:,:]
	#Positions
	X, Y, Z = size(ρ)
	FOVx = (X-1)*Δx[1] #[m]
	FOVy = (Y-1)*Δx[2] #[m]
	FOVz = (Z-1)*Δx[3] #[m]
	xx = [(-FOVx/2:Δx[1]:FOVx/2)...;]
	yy = [(-FOVy/2:Δx[2]:FOVy/2)...;;]
	zz = [(-FOVz/2:Δx[3]:FOVz/2)...;;;]
	x = xx*1 .+ yy*0 .+ zz*0 .+ offset[1]	#spin x coordinates
	y = xx*0 .+ yy*1 .+ zz*0 .+ offset[2]	#spin y coordinates
	z = xx*0 .+ yy*0 .+ zz*1 .+ offset[3]	#spin z coordinates
	v = 0 # m/s

	phantom = Phantom(name=basename(filename),
					x=    x[mask],
					y=    y[mask],
					z=    z[mask],
					ρ=    ρ[mask],
					T1=	 T1[mask],
					T2=	 T2[mask],
					T2s=T2s[mask],
					Δw=	 Δw[mask],
					ux=(x,y,z,t)->v*t
					)
	phantom
end
