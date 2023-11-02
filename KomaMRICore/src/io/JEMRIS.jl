"""
    obj = read_phantom_jemris(filename)

Returns the Phantom struct from a JEMRIS phantom file `.h5`.

# Arguments
- `filename`: (`::String`) the absolute or relative path of the phantom file `.h5`

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj_file = joinpath(dirname(pathof(KomaMRI)), "../examples/2.phantoms/brain.h5")

julia> obj = read_phantom_jemris(obj_file)

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
	xx = reshape((-FOVx/2:Δx[1]:FOVx/2),:,1,1) #[(-FOVx/2:Δx[1]:FOVx/2)...;]
	yy = reshape((-FOVy/2:Δx[2]:FOVy/2),1,:,1) #[(-FOVy/2:Δx[2]:FOVy/2)...;;]
	zz = reshape((-FOVz/2:Δx[3]:FOVz/2),1,1,:) #[(-FOVz/2:Δx[3]:FOVz/2)...;;;]
	x = xx*1 .+ yy*0 .+ zz*0 .+ offset[1]	#spin x coordinates
	y = xx*0 .+ yy*1 .+ zz*0 .+ offset[2]	#spin y coordinates
	z = xx*0 .+ yy*0 .+ zz*1 .+ offset[3]	#spin z coordinates
	v = 0 # m/s

	obj = Phantom{Float64}(
        name = basename(filename),
		x = x[mask],
		y = y[mask],
		z = z[mask],
		ρ = ρ[mask],
		T1 = T1[mask],
		T2 = T2[mask],
		T2s = T2s[mask],
		Δw = Δw[mask],
		ux = (x,y,z,t)->v*t,
	)
	return obj
end
