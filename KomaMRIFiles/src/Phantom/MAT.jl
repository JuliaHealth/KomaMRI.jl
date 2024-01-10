"""
phantom = read_phantom_MAT(filename)

Reads a (.mat) file and creates a Phantom structure from it
"""
function read_phantom_MAT(folder::String; ss::Int=1, Δx=1)
	ρ  = matread(folder*"mapaPD.mat")["PD"]
	T1 = matread(folder*"mapaT1.mat")["T1"]
	T2 = matread(folder*"mapaT2.mat")["T2"]

	# Subsampling
	ρ  = ρ[1:ss:end,1:ss:end,1:ss:end]
	T1 = T1[1:ss:end,1:ss:end,1:ss:end]
	T2 = T2[1:ss:end,1:ss:end,1:ss:end]

	# Clip outliers
	T1_percentile = KomaMRIBase.percentile(T1[:],99)
	T1[T1.>=T1_percentile].= T1_percentile;

	T2_percentile = KomaMRIBase.percentile(T2[:],99)
	T2[T2.>=T2_percentile].= T2_percentile;

	ρ_percentile = KomaMRIBase.percentile(ρ[:],99)
	ρ[ρ.>=ρ_percentile].= ρ_percentile;

	# Normalize ρ between 0 and 1
	mini, maxi = extrema(ρ)
    ρ = (ρ .- mini) ./ (maxi - mini);

	# Take small ρ values to 0
	thresh = 0.2
	ρ[ρ.<=thresh] .= 0

	# Convert miliseconds into seconds
	T1 .*= 1e-3
	T2 .*= 1e-3

	Δx *= 1e-3*ss  

	M, N, L = size(ρ)

	FOVx = (M-1)*Δx
	FOVy = (N-1)*Δx
	FOVz = (L-1)*Δx

	println("Phantom dimensions: ($FOVx x $FOVy x $FOVz) m")

	xx = reshape(-FOVx/2:Δx:FOVx/2,M,1,1) 
    yy = reshape(-FOVy/2:Δx:FOVy/2,1,N,1) 
	zz = reshape(-FOVz/2:Δx:FOVz/2,1,1,L) 

	# Grid
    # X = 1*xx .+ 0*yy .+ 0*zz
	# Y = 0*xx .+ 1*yy .+ 0*zz
	# Z = 0*xx .+ 0*yy .+ 1*zz

	X = matread(folder*"X.mat")["X"] .* 1e-3
	Y = matread(folder*"Y.mat")["Y"] .* 1e-3
	Z = matread(folder*"Z.mat")["Z"] .* 1e-3

	X  = X[1:ss:end,1:ss:end,1:ss:end]
	Y  = Y[1:ss:end,1:ss:end,1:ss:end]
	Z  = Z[1:ss:end,1:ss:end,1:ss:end]



	println("Number of spins: ", length(ρ[ρ.!=0]))

	name = split(folder,"/")[end]

	deltaX = matread(folder*"deltaX.mat")["deltaX"] .* 1e-3
	deltaY = matread(folder*"deltaY.mat")["deltaY"] .* 1e-3
	deltaZ = matread(folder*"deltaZ.mat")["deltaZ"] .* 1e-3

	deltaX  = deltaX[1:ss:end, 1:ss:end, 1:ss:end, :]
	deltaY  = deltaY[1:ss:end, 1:ss:end, 1:ss:end, :]
	deltaZ  = deltaZ[1:ss:end, 1:ss:end, 1:ss:end, :]

	K = size(deltaX)[4] + 1

	Δx = zeros(length(ρ[ρ.!=0]),K-1)
	Δy = zeros(length(ρ[ρ.!=0]),K-1)
	Δz = zeros(length(ρ[ρ.!=0]),K-1)

	for i in 1:K-1
		Δx[:,i] = deltaX[:,:,:,i][ρ.!=0]
		Δy[:,i] = deltaY[:,:,:,i][ρ.!=0]
		Δz[:,i] = deltaZ[:,:,:,i][ρ.!=0]
	end


	resetmag = BitMatrix(zeros(length(ρ[ρ.!=0]),K))

	phantom = Phantom{Float64}(name="Heart XCAT",
                          		x=X[ρ.!=0],
                          		y=Y[ρ.!=0],
                          		z=Z[ρ.!=0],
                          		ρ=ρ[ρ.!=0],
                          		T1=T1[ρ.!=0],
                          		T2=T2[ρ.!=0],
								motion = ArbitraryMotion(	[1.0],
															K,
															Δx,
															Δy,
															Δz,
															resetmag))
	phantom
end