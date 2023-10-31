"""
Simulate cardiac cine adquisition and reconstruction. 
Output are cine frames from different cardiac phases
"""
cardiac_cine(FOV::Float64,heart_rate::Int,N_phases::Int,N::Int,obj::Phantom,sys::Scanner) = begin
	RR = 60/heart_rate
	TR = RR/N_phases
	α = 15

	# Sequence
	seq = Sequence()
	base_seq =  PulseDesigner.bSSFP(FOV, N, TR, α, sys)
	for i in 1:N
		for j in 1:N_phases
			seq += base_seq[6*(i-1).+(1:6)]
		end
	end
	seq = seq[2:end]

	# Simulation
	raw_signal = simulate(obj, seq, sys)

	# Reconstruction
	recParams = Dict{Symbol,Any}(:reco=>"direct")
	Nx, Ny = raw_signal.params["reconSize"][1:2]
	recParams[:reconSize] = (Nx, Ny)
	recParams[:densityWeighting] = false

	acqData = AcquisitionData(raw_signal)

	frames = []
	for i in 1:N_phases
		acqAux = copy(acqData)
		range = reduce(vcat,[j*(N*N_phases).+((i-1)*N.+(1:N)) for j in 0:N-1])

		# Kdata
		acqAux.kdata[1] = reshape(acqAux.kdata[1][range],(N^2,1))

		# Traj
		acqAux.traj[1].circular = false

		acqAux.traj[1].nodes = acqAux.traj[1].nodes[:,range]
		acqAux.traj[1].nodes = acqAux.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqAux.traj[1].nodes[:]))

		acqAux.traj[1].numProfiles = N
		acqAux.traj[1].times = acqAux.traj[1].times[range]

		# subsampleIndices
		acqAux.subsampleIndices[1] = acqAux.subsampleIndices[1][1:N^2]

		# Reconstruction
		aux = @timed reconstruction(acqAux, recParams)
		image  = reshape(aux.value.data,Nx,Ny,:)
		image_aux = abs.(image[:,:,1])
		# image_aux = round.(UInt8,255*(image_aux./maximum(image_aux)))

		push!(frames,image_aux)
	end
	frames
end