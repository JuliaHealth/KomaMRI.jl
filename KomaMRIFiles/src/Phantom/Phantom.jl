"""
phantom = read_phantom(filename)

Reads a (.phantom) file and creates a Phantom structure from it
"""
function read_phantom(filename::String) 
	# ----------------------------------------------
	function read_param(param::HDF5.Group)

		if "type" in HDF5.keys(attrs(param))
			type = attrs(param)["type"]
		
			if     type == "Explicit"
				values = read(param["values"])
			elseif type == "Indexed"
				index = read(param["values"]) 
				if Ns == length(index)
					table = read(param["table"])
					N = read_attribute(param,"N")
					if N == length(table)
						values = table[index]
					else
						print("Error: $(label) table dimensions mismatch")
					end
				else
					print("Error: $(label) vector dimensions mismatch")
				end
			elseif type == "Default"
				values = "Default"
			end
		else
			values = read(param["values"])
		end
	
		values
	end
	# ----------------------------------------------

	fid = HDF5.h5open(filename,"r")

	name    = read_attribute(fid,"Name")
	version = read_attribute(fid,"Version")
    dims    = read_attribute(fid,"Dims")
    dynamic = Bool(read_attribute(fid,"Dynamic"))
    Ns      = read_attribute(fid,"Ns")

	obj = Phantom(name=name,
				  x=zeros(Ns))

	# Position and contrast
	for key in ["position","contrast"]
		group = fid[key]
		for label in HDF5.keys(group)
			param = group[label]
			values = read_param(param)
			if values != "Default"
				setfield!(obj,Symbol(label),values)
			end
		end
	end

	# Motion
	if dynamic
		motion = fid["motion"]
		model = read_attribute(motion,"model")
		if model == "Simple"
			# name = motion["name"]
			# ux, uy, uz = get_simple_motion(name)
			
			# SimpleMotion(
			# 	ux = ux,
			# 	uy = uy,
			# 	uz = uz
			# )
		elseif model == "Arbitrary"
			segments = motion["segments"]
			N = read_attribute(segments, "N")
			K = read_attribute(segments, "K")
			dur = read(segments["dur"])

			Δx = zeros(Ns,K-1)
			Δy = zeros(Ns,K-1)
			Δz = zeros(Ns,K-1)

			resetmag = zeros(Ns,K)

			for key in HDF5.keys(motion)
				if key != "segments"
					values = read_param(motion[key])
					if 	   key == "motionx"
						Δx = values
					elseif key == "motiony"
						Δy = values
					elseif key == "motionz"
						Δz = values
					elseif key == "resetmag"
						resetmag = Bool.(values)
					end
				end
			end

			obj.motion = ArbitraryMotion(	
							Float64.(dur),
							K,
							Float64.(Δx),
							Float64.(Δy),
							Float64.(Δz),
							resetmag
						)
		end
	end

	close(fid)
	return obj
end


"""
write_phantom(ph,filename)

Writes a (.phantom) file from a Phantom struct.
"""
# By the moment, only "Explicit" type 
# is considered when writing .phantom files
function write_phantom(obj::Phantom,filename::String) 
	# Create HDF5 phantom file
	fid = h5open(filename,"w")

	# Root attributes
	HDF5.attributes(fid)["Version"] = "1.0"
	HDF5.attributes(fid)["Name"] = obj.name
	HDF5.attributes(fid)["Ns"] = length(obj.x)
	dims = get_dims(obj)
	HDF5.attributes(fid)["Dims"] = sum(dims)
	dynamic = is_dynamic(obj.motion)
	HDF5.attributes(fid)["Dynamic"] = Int(dynamic)     # 0=False, 1=True

	fields = fieldnames(Phantom)[2:end]

	# Spin initial positions
	pos = create_group(fid,"position")
	for i in 1:3
		if Bool(dims[i])
			create_group(pos,String(fields[i]))["values"] = getfield(obj,fields[i])
		end
	end

	# Contrast (Rho, T1, T2, T2s Deltaw)
	contrast = create_group(fid,"contrast")
	for i in 4:8
		param = create_group(contrast,String(fields[i]))
		HDF5.attributes(param)["type"] = "Explicit"
		param["values"] = getfield(obj,fields[i])
	end

	# Motion
	if dynamic
		motion = create_group(fid,"motion")
		if typeof(obj.motion) <: SimpleMotion
			# HDF5.attributes(motion)["model"] = "Simple"
			# tmp_path = tempname()*".jld2"
			# JLD2.save(tmp_path,Dict("ux" => obj.motion.ux,
			# 						"uy" => obj.motion.uy,
			# 						"uz" => obj.motion.uz))
			# tmp = h5open(tmp_path,"r")
			# for key in keys(tmp)
			# 	copy_object(tmp,key,motion,key)
			# end

		elseif typeof(obj.motion) <: ArbitraryMotion
			HDF5.attributes(motion)["model"] = "Arbitrary"

			segments = create_group(motion, "segments")
			HDF5.attributes(segments)["N"] = length(obj.motion.dur) 
			HDF5.attributes(segments)["K"] = obj.motion.K
			segments["dur"] = obj.motion.dur

			itp = get_itp_functions(obj.motion)[1]
			is_mov_on = vcat((itp .!== nothing),is_fluid(obj.motion)) 
			mov_dims = ["motionx","motiony","motionz","resetmag"] 
			Δ = fieldnames(ArbitraryMotion)[3:6]
			for i in 1:4
				if is_mov_on[i]
					motion_i = create_group(motion,mov_dims[i])
					HDF5.attributes(motion_i)["type"] = "Explicit"
					values = i != 4 ? getfield(obj.motion,Δ[i]) : Int.(getfield(obj.motion,Δ[i]))
					motion_i["values"] = values
				end
			end
		end
	end

	close(fid)
end
