"""
    phantom = read_phantom_MRiLab(filename)

Returns the Phantom struct from a MRiLab phantom file `.mat`.

# Arguments
- `filename`: (`::String`) the absolute or relative path of the phantom file `.mat`

# Returns
- `phantom`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj_file = joinpath(dirname(pathof(KomaMRI)), "../examples/2.phantoms/brain.mat")

julia> obj = read_phantom_MRiLab(obj_file)

julia> plot_phantom_map(obj, :ρ)
```
"""
function read_phantom_MRiLab(filename; B0=1.5, offset=[0,0,0], FRange_filename="")
    data = MAT.matread(filename)["VObj"]

	Δx = [data["XDimRes"], data["YDimRes"], data["ZDimRes"]] #[m]
	mask = data["Rho"] .!= 0
	if FRange_filename!=""
		data_mask = MAT.matread(FRange_filename)["VMag"]
		mask .*= data_mask["FRange"]
	end
	#Maps
	ρ =   data["Rho"]
	T1 =  data["T1"]
	T2 =  data["T2"]
	T2s = data["T2Star"]
    CS = data["ChemShift"] #[Hz/T]
	Δw = 2π .* CS .* B0 .* ones(size(ρ))
	#Positions
	X, Y, Z = size(ρ)
	FOVx = (X-1)*Δx[1] #[m]
	FOVy = (Y-1)*Δx[2] #[m]
	FOVz = (Z-1)*Δx[3] #[m]
    # [row,col,layer]=size(VOex.Mz);
    # VVar.ObjLoc =
            # [((col+1)/2)*VOex.XDimRes;
            # ((row+1)/2)*VOex.YDimRes ;
            # ((layer+1)/2)*VOex.ZDimRes]; % Set matrix center as Object position for motion simulation
    # VVar.ObjTurnLoc =
            # [((col+1)/2)*VOex.XDimRes;
            # ((row+1)/2)*VOex.YDimRes ;
            # ((layer+1)/2)*VOex.ZDimRes]; % Set matrix center as Object origin for motion simulation
    xx = reshape((-FOVx/2:Δx[1]:FOVx/2),:,1,1) #[(-FOVx/2:Δx[1]:FOVx/2)...;]
	yy = reshape((-FOVy/2:Δx[2]:FOVy/2),1,:,1) #[(-FOVy/2:Δx[2]:FOVy/2)...;;]
	zz = reshape((-FOVz/2:Δx[3]:FOVz/2),1,1,:) #[(-FOVz/2:Δx[3]:FOVz/2)...;;;]
	x = xx*1 .+ yy*0 .+ zz*0 .+ offset[1]	#spin x coordinates
	y = xx*0 .+ yy*1 .+ zz*0 .+ offset[2]	#spin y coordinates
	z = xx*0 .+ yy*0 .+ zz*1 .+ offset[3]	#spin z coordinates
	v = 0 # m/s

	phantom = Phantom(name=basename(filename),
					x=    y[mask],
					y=   -x[mask],
					z=   -z[mask],
					ρ=    ρ[mask],
					T1=	 T1[mask],
					T2=	 T2[mask],
					T2s=T2s[mask],
					Δw=	 Δw[mask],
					ux=(x,y,z,t)->v*t
					)
	phantom
end
