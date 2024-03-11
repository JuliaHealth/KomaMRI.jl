"""
Simulate cardiac cine adquisition and reconstruction. 
Output are cine frames from different cardiac phases
"""
cardiac_cine(FOV::Float64,heart_rate::Int,N_phases::Int,N::Int,obj::Phantom,sys::Scanner;
			 Δf=0,
			 flip_angle = 10,
			 TR = (60/heart_rate)/N_phases,
			 dummy_cycles = 0,
			 tagging::Bool=false) = begin

	RR = 60/heart_rate
	α = flip_angle

	# Sequence
	global seq = Sequence()


	# Tagging ----------------------------------
	tag = Sequence()
	if tagging

		# SPAMM
		hard_flip(T,α,sys) = begin
			B1 = α/(360*γ*T) 
			return PulseDesigner.RF_hard(B1, T, sys)
		end

		T_RF = 0.5e-3

		EX_22 = hard_flip(T_RF,22.5,sys)
		EX_45 = hard_flip(T_RF,45,sys)
		EX_90 = hard_flip(T_RF,90,sys)

		A = 3.5e-3
		T = 0.7e-3
		ζ = A / sys.Smax

		GR_x = Sequence(reshape([Grad(A,T,ζ);
								Grad(0,0);
								Grad(0,0)],(3,1)))

		GR_y = Sequence(reshape([Grad(0,0);
								Grad(A,T,ζ);
								Grad(0,0)],(3,1)))

		spamm_x =   EX_45 +
					GR_x +
					EX_45 +
					5*GR_x

		spamm_y =   EX_45 +
					GR_y +
					EX_45 +
					5*GR_y

		tag = spamm_x + spamm_y + Delay(10e-3)


		# DANTE 
		"""
		dante = PulseDesigner.RF_train(8, 1e-4, 1e-3, 25, sys; G = [0,4e-3,0])

		tag = dante + rotz(π/2)*dante + Delay(1e-3)
		"""
	end
	# -------------------------------------------------

	prospective = true
	if TR == RR/N_phases
		TR -= sum(dur(tag))/N_phases
		prospective = false
	end

	base_seq =  PulseDesigner.bSSFP(FOV, N, TR, α, sys; Δf=Δf)


	for i in 0:N-1
		seq += tag

		line = base_seq[6*i .+ (1:6)]
		dummy_cycle = copy(line)
		dummy_cycle.ADC = [ADC(0,0) for i in 1:length(dummy_cycle)]
		# dummy_cycle.GR[1,:] = [Grad(0,0) for i in 1:length(dummy_cycle)]
		# dummy_cycle.GR[2,:] = [Grad(0,0) for i in 1:length(dummy_cycle)]

		for j in 1:dummy_cycles
			seq += dummy_cycle
		end
		for j in 1:N_phases
			seq += line
		end
		if prospective
			dead_space = RR - (N_phases*sum(dur(base_seq[6*i .+ (1:6)])) + sum(dur(tag)))
			seq += Delay(dead_space)
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