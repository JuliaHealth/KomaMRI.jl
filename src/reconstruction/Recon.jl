ifftc(x;dims=[1,2])=fftshift(ifft(ifftshift(x,dims),dims),dims)*prod(size(x)[dims])
fftc(x;dims=[1,2]) =fftshift(fft(ifftshift(x,dims),dims),dims)

#DCF 
"""
Hoge, R.D., Kwan, R.K.S. and Bruce Pike, G. (1997), Density compensation functions for spiral MRI. Magn. Reson. Med., 38: 117-128. https://doi.org/10.1002/mrm.1910380117"""
function dcf_nufft(seq)
    _, k = get_kspace(seq)
	t_adc = get_sample_times(seq) |> t-> convert(Vector{Float64},t)
	Gx, Gy, Gz = get_grads(seq, t_adc)
	dcf = abs.(Gx .* k[:,1] .+ Gy .* k[:,2] .+ Gz .* k[:,3]) # dcf = γ|<k(t), g(t)>|
	dcf ./ maximum(dcf[:])
end

# using MRIReco: AcquisitionData
# function reconstruction(signal, recParams)
#     #Param unpack
#     Nx =  get(recParams, "Nx",  101)				|> x->floor(Int,x)
# 	Ny =  get(recParams, "Ny",  recParams["Nx"]) 	|> x->floor(Int,x)
# 	Nz =  get(recParams, "Nz",  1) 					|> x->floor(Int,x)
# 	name = get(recParams, "Name", "")
#     recon = get(recParams, "recon", "MIRT")

# 	println("")
# 	@info "Running reconstruction..."
# 	@time begin
# 	#K-data, only 2D for now
# 	signal_phase = signal .* get_sample_phase_compensation(seq) #ADC phase compensation
# 	kdata = reshape(signal_phase, Nx,Ny,Nz)
# 	image = zeros(2,2,1)
# 	kspace = zeros(2,2,1)
# 	if recon != "skip" 
# 		cartesian = !(occursin("rad", name) || occursin("spiral", name))
# 		if recon == "MRIReco"
# 			#MRIReco.jl, WIP
# 			#Getting kspace coordinates
# 			times = KomaMRI.get_sample_times(seq)
# 			_, traj = KomaMRI.get_kspace(seq)
# 			TE = 0
# 			AQ = maximum(seq[KomaMRI.is_ADC_on.(seq)].ADC.T)
# 			tr = Trajectory(name, traj', times, TE, AQ, Ny, Nx, Nz, cartesian, !cartesian)
# 			#RawData
# 			sequenceInfo = Dict([Symbol(key) => value for (key,value) = seq.DEF]...)
# 			kdata = [zeros(ComplexF64,size(nodes,2),nc) for echo=1:1, slice=1:nz, rep=1:1]
# 			raw = AcquisitionData(sequenceInfo,tr,kdata)
# 			#(1. dim k-spacenodes, 2. dim coils) the outer dims describe: 
# 			#1. dim echoes, 2. dim slices, 3. dim repetitions
# 			println(raw)
# 			# sequenceInfo::Dict{Symbol,Any}
# 			# traj::Vector{Trajectory}
# 			# kdata::Array{Matrix{ComplexF64},3}
# 			# subsampleIndices::Vector{Vector{Int64}}
# 			# encodingSize::Vector{Int64}
# 			# fov::Vector{Float64}
# 		end
# 		if recon=="MIRT"
# 			#MIRT.jl
# 			#Getting kspace coordinates
# 			_, traj = get_kspace(seq)
# 			k2d = k[:,1:2]
# 			#Normalize
# 			kmax = maximum(abs.(k2d[:]))
# 			ω = π * k2d / kmax
# 			#Number of samples and Operator definition, TODO: consider vector FOV vectors
# 			if non_cartesian
# 				N = 2 * round(Int, kmax * seq.DEF["FOV"] )
# 				Nu = Nv = N
# 				dcf = dcf_nufft(seq)
# 			else
# 				Nu = Nx
# 				Nv = Ny
# 				dcf = 1
# 			end
# 			image = zeros(ComplexF64,Nu,Nv,Nz)
# 			#MIRT.jl NUFFT
# 			for i = 1:Nz
# 				n_shift=round.(Int, [Nu/2; Nv/2])
# 				p = nufft_init(ω, (Nu,Nv); n_shift, pi_error=false)
# 				#Reconstruction
# 				kdata_i = kdata[:,:,i]
# 				print(size(kdata_i))
# 				image[:,:,i] = p.A' * (dcf .* kdata_i[:])
# 			end
# 			println(size(image))
# 			kspace = fftc(image, dims=[1,2])
# 		end
# 	end
# 	end
# 	(image, kspace)
# end