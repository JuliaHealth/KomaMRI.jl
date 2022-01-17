ifftc(x::Array{Complex{Float64},2})=fftshift(ifft(ifftshift(x)))*prod(size(x))
fftc(x::Array{Complex{Float64},2})=fftshift(fft(ifftshift(x)))
function reconstruction(signal, recParams)
    #Param unpack
    Nx =  get(recParams, "Nx",  101)
	Ny =  get(recParams, "Ny",  recParams["Nx"])
	epi = get(recParams, "EPI", false)
    recon = get(recParams, "recon", "fft")
	#K-data, only 2D for now
	kdata = reshape(copy(signal),(Nx,Ny)) #Turning into kspace image
	if epi #Flip in freq-dir for the EPI
		kdata[:,2:2:Ny] = kdata[Nx:-1:1,2:2:Ny] 
	end 
	kdata = convert(Array{Complex{Float64},2},kdata)
	if recon != "skip" 
		#Recon, will be replaced to call MRIReco.jl
		if recon == "fft"
			image = ifftc(kdata)
		end
        image
	else
		kspace = log.(abs.(kdata) .+ 1)
        kspace
	end
end
