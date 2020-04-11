ifftc(x::Array{Complex{Float64},2})=fftshift(ifft(ifftshift(x)))*prod(size(x))
fftc(x::Array{Complex{Float64},2})=fftshift(fft(ifftshift(x)))
