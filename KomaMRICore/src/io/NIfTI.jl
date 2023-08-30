"""
phantom = read_phantom_NIfTI(folder; ss)
Reads a (.nii) file and creates a Phantom structure from it
"""

function read_phantom_NIfTI(folder; ss=1::Int)
	T1_path = folder*"/T1map.nii.gz";
	T2_path = folder*"/T2map.nii.gz";
	PD_path = folder*"/PDmap.nii.gz";

	T1_ni = niread(T1_path);
	T2_ni = niread(T2_path);
	PD_ni = niread(PD_path);

	T1_data = T1_ni.raw
	T2_data = T2_ni.raw
	ρ = PD_ni.raw

	M, N, L = size(T1_data)

	# Subsampling
	T1_data = T1_data[1:ss:end,1:ss:end,1:ss:end]
	T2_data = T2_data[1:ss:end,1:ss:end,1:ss:end]
	ρ = ρ[1:ss:end,1:ss:end,1:ss:end]

	# Clip outliers
	T1_percentile = percentile(T1_data[:],99)
	T1_data[T1_data.>=T1_percentile].= T1_percentile;

	T2_percentile = percentile(T2_data[:],99)
	T2_data[T2_data.>=T2_percentile].= T2_percentile;

	ρ_percentile = percentile(ρ[:],99)
	ρ[ρ.>=ρ_percentile].= ρ_percentile;

	# Normalize ρ between 0 and 1
	mini, maxi = extrema(ρ)
    ρ = (ρ .- mini) ./ (maxi - mini);

	# Take small ρ values to 0
	thresh = 0.15
	ρ[ρ.<=thresh] .= 0

	Δx = 1e-3*ss  # Each voxel is 1mm x 1mm x 1mm.
				  # Voxels have the same size in the three dimensions (cubic voxels)

	M, N, L = size(T1_data)

	FOVx = (M-1)*Δx
	FOVy = (N-1)*Δx
	FOVz = (L-1)*Δx

	println("Phantom dimensions: ($FOVx x $FOVy x $FOVz) m")

	xx = reshape(-FOVx/2:Δx:FOVx/2,M,1,1) #spin coordinates
    yy = reshape(-FOVy/2:Δx:FOVy/2,1,N,1) #spin coordinates
	zz = reshape(-FOVz/2:Δx:FOVz/2,1,1,L) #spin coordinates
	# Grid
    X = 1*xx .+ 0*yy .+ 0*zz
	Y = 0*xx .+ 1*yy .+ 0*zz
	Z = 0*xx .+ 0*yy .+ 1*zz

	Δw = zeros(M,N,L)

	# Convert miliseconds into seconds
	T1_data .*= 1e-3
	T2_data .*= 1e-3

	println("Number of spins: ", length(X[ρ.!=0]))

	name = split(folder,"/")[end]

	phantom = Phantom{Float64}(name=name,
                          		x=X[ρ.!=0],
                          		y=Y[ρ.!=0],
                          		z=Z[ρ.!=0],
                          		ρ=ρ[ρ.!=0],
                          		T1=T1_data[ρ.!=0],
                          		T2=T2_data[ρ.!=0],
                          		Δw=Δw[ρ.!=0]
                          		)
	phantom
end