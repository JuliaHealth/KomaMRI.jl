ifftc(x::Array{Complex{Float64},2})=fftshift(ifft(ifftshift(x)))*prod(size(x))
fftc(x::Array{Complex{Float64},2})=fftshift(fft(ifftshift(x)))

# kdata = convert(Array{Complex{Float64},2},transpose(kdata_raw))
# Ploting recon
# rec = ifftc(kdata)
# plot_recon(kdata,rec,Δx_pix,Δx_pix,"./Results/SIMPLE_asym_recon_t"*string(i),title)
