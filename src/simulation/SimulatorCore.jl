##########################
## SIMULATION FUNCTIONS ##
##########################

#GPU related functions
print_gpus() = begin
	println( "$(length(devices())) CUDA capable device(s)." )
	for (i,d) = enumerate(devices())
		u = i == 1 ? "*" : " "
		println( "  ($(i-1)$u) $(name(d))")
	end
end

"""
Implementation in multiple threads. Separating the spins in N_parts.
"""
function run_spin_precession_parallel(obj::Phantom,seq::Sequence, t::Array{Float64,1}, Δt::Array{Float64,1};
	M0::Array{Mag,1}, Nthreads::Int=Hwloc.num_physical_cores(), gpu::Function)

	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	#Put times as row vector
	t = reshape(t,1,Nt)
	Δt = reshape(Δt,1,NΔt)

	S = zeros(ComplexF64, Nt)
	
	parts = kfoldperm(Ns, Nthreads, type="ordered") 

	@threads for p ∈ parts
		@inbounds aux, M0[p] = run_spin_precession(obj[p],seq,t,Δt; M0=M0[p], gpu)
		S .+= aux
		aux = nothing
	end
    S, M0
end

"""
	run_spin_precession(obj,seq,t)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`.
It calculates S(t) = ∫ ρ(x,t) exp(- t/T2(x,t) ) exp(- 𝒊 ϕ(x,t)) dx.
"""
function run_spin_precession(obj::Phantom, seq::Sequence, t::Array{Float64,2}, Δt::Array{Float64,2};
	M0::Array{Mag,1}, gpu::Function)

	𝒊 = 1im; Random.seed!(1)
	T = sum(Δt) #Total length, used for signal relaxation
	t = [t t[end]+Δt[end]] #For trapezoidal integration <---------------------------
    sz = size(obj)
	Nsz = length(sz)

	#DIFFUSION, disabled while I think how to do it efficiently
	# if any(obj.Dλ1 .+ obj.Dλ2 .!= 0)  #No diff optimization
	# 	#TODO: I need to add diff displacement η story between blocks (like ϕ0) <-- This is already taken care of in M0
	# 	η1 = randn(sz...,Nt) |> gpu
	# 	η2 = randn(sz...,Nt) |> gpu
	# 	Dθ = obj.Dθ |> gpu
	# 	Dλ1 = obj.Dλ1 |> gpu
	# 	Dλ2 = obj.Dλ2 |> gpu
	# 	ηx = sqrt.(2Δt .* Dλ1) .* η1 |> gpu
	# 	ηy = sqrt.(2Δt .* Dλ2) .* η2 |> gpu
	# 	ηxp = cumsum(ηx .* cos.(Dθ) .- ηy.*sin.(Dθ), dims=Nsz+1) |> gpu
	# 	ηyp = cumsum(ηy .* cos.(Dθ) .+ ηx.*sin.(Dθ), dims=Nsz+1) |> gpu
	# 	ηzp = 0
	# else
		ηxp = 0
		ηyp = 0
		ηzp = 0
	# end

	#SCANNER
    Gx, Gy, Gz = get_grads(seq,t)
	Gx = Gx |> gpu
	Gy = Gy |> gpu 
	Gz = Gz |> gpu 
	#SIMULATION
	Mxy = abs.(M0)	|> gpu
	ϕ0 = angle.(M0)	|> gpu
	x0 = obj.x		|> gpu
	y0 = obj.y		|> gpu
	z0 = obj.z		|> gpu
	tp = t[:,1:end-1] .- t[1]	|> gpu # t' = t - t0
	t = t			|> gpu
	Δt = Δt			|> gpu
    xt = x0 .+ obj.ux(x0,y0,z0,t) .+ ηxp |> gpu
	yt = y0 .+ obj.uy(x0,y0,z0,t) .+ ηyp |> gpu
	zt = z0 .+ obj.uz(x0,y0,z0,t) .+ ηzp |> gpu
	#ACQ OPTIMIZATION
    if is_ADC_on(seq, Array(t)) 
		ϕ =  ϕ0 .- (2π*γ) .* cumtrapz(Δt, xt.*Gx .+ yt.*Gy .+ zt.*Gz )
	else
		ϕ =  ϕ0 .- (2π*γ) .* trapz(Δt, xt.*Gx .+ yt.*Gy .+ zt.*Gz )
	end
	#Mxy preccesion and relaxation
	Δw = obj.Δw  |> gpu #Need to add a component here to model scanner's dB0(xt,yt,zt)
	T2 = obj.T2  |> gpu
	Mxy =  Mxy .* exp.(𝒊.*(ϕ .- Δw.*tp) .- tp./T2 ) #This assumes Δw constant in time
	#ACQUIRED SIGNAL
	S = sum(Mxy, dims=1:Nsz)[:] #<--- TODO: add coil sensitivities
	#Mz relaxation
	T1 = obj.T1		|> gpu
	Mz0 = obj.ρ		|> gpu
	Mz = M0.z		|> gpu
	Mz =  Mz .* exp.(-T./T1) .+ Mz0 .* ( 1 .- exp.(-T./T1) )
	# END
	Mxy = Array(Mxy)
	Mz = Array(Mz)
	M0 = Mag.(Mxy[:,end], Mz) #Saving the last magnetization
    Array(S), M0 #Singal, M0
end

run_spin_excitation_parallel(obj, seq, t::Array{Float64,1}, Δt::Array{Float64,1}; 
	M0::Array{Mag,1}, Nthreads::Int=Hwloc.num_physical_cores(), gpu::Function) = begin
	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	#Put times as row vector
	t = reshape(t,1,Nt)
	Δt = reshape(Δt,1,NΔt)
	
	parts = kfoldperm(Ns, Nthreads, type="ordered") 

	@threads for p ∈ parts
		@inbounds M0[p] = run_spin_excitation(obj[p],seq,t,Δt; M0=M0[p], gpu)
	end
    M0
end

run_spin_excitation(obj, seq, t::Array{Float64,2}, Δt::Array{Float64,2}; 
	M0::Array{Mag,1}, gpu::Function) = begin
	#SCANNER
	B1, Δf_rf  = get_rfs(seq,t)
    Gx, Gy, Gz = get_grads(seq,t)
	B1 = B1 |> gpu
	Gx = Gx |> gpu
	Gy = Gy |> gpu 
	Gz = Gz |> gpu 
	#SIMULATION
	x0 = obj.x		|> gpu
	y0 = obj.y		|> gpu
	z0 = obj.z		|> gpu
	t = t			|> gpu
	Δt = Δt			|> gpu
    xt = x0 .+ obj.ux(x0,y0,z0,t)		|> gpu
	yt = y0 .+ obj.uy(x0,y0,z0,t)		|> gpu
	zt = z0 .+ obj.uz(x0,y0,z0,t)		|> gpu
	ΔB0 = obj.Δw./(2π*γ) .- Δf_rf./γ	|> gpu # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
	Bz = (Gx.*xt .+ Gy.*yt .+ Gz.*zt) .+ ΔB0	#<-- TODO: This line is very slow, FIX!?
	B = sqrt.(abs.(B1).^2. .+ abs.(Bz).^2.)		
	φ = -2π*γ * (B .* Δt) # angle of rotation 
	B[B.==0] .= 1e-17; # removes problems when dividing by φ
	Qt = Q.(Array(φ), Array(B1./B), Array(Bz./B)) #TODO: remove the gpu->array step
	Qf = prod( Qt , dims=2 )[:] # equivalent rotation
	M0 = Qf .* M0 #TODO: This and the relaxation are computed in the CPU for now
	#Relaxation
	T = sum(Δt) #Total length, used for signal relaxation
	for (m, ρ, T1, T2) = zip(M0, obj.ρ, obj.T1, obj.T2)
		m.xy *= @fastmath exp(-T/T2)
		m.z = @fastmath m.z*exp(-T/T1) + ρ*(1-exp(-T/T1))
	end
 	M0
end

"""Divides time steps in N_parts blocks. Decreases RAM usage in long sequences."""
function run_sim_time_iter(obj::Phantom,seq::Sequence, t::Array{Float64,1}, Δt; 
								Nblocks::Int=16, Nthreads::Int=Hwloc.num_physical_cores(), gpu::Function, w=nothing)
	Nt, Ns = length(t), prod(size(obj))
	blink_window = w !== nothing
	S = zeros(ComplexF64, Nt) #Only one coil for now, TODO: change for more coils to Nt x Coils
	M0 = Mag(obj,:z)
	breaks = get_breaks_in_RF_key_points(seq,t)
    parts = kfoldperm(Nt,Nblocks;type="ordered",breaks)
	Nblocks = length(parts)
	println("Dividing simulation in Nblocks=$Nblocks")
	println("Starting simulation with Nspins=$Ns and Nt=$Nt")
	#Perturbation of spins' position to reduce spurious echoes (?)
	#Test convert T2* to ΔBz with Lorentzian distribution (?)
	# R2prime = 1 ./obj.T2s .- 1 ./ obj.T2 #1/T2* = 1/T2 + 1/T2' and 1/T2' = γΔB
	#obj_p.Δw .+= something
	#TODO: transform suceptibility χ to Δω, for each time-block with FMM-like technique O(nlogn).
	rfs = 0
	pp = Progress(Nblocks)
	for (block, p) = enumerate(parts)
		if is_RF_on(seq, t[p]) && !is_ADC_on(seq, t[p]) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
			@inbounds M0  = run_spin_excitation_parallel(obj, seq, t[p], Δt[p]; M0, Nthreads, gpu)
			rfs += 1
		else
			@inbounds S[p], M0 = run_spin_precession_parallel(obj, seq, t[p], Δt[p]; M0, Nthreads, gpu)
		end
		#Update progress
		next!(pp, showvalues = [(:simulated_blocks, block), (:rf_blocks,rfs)])
		if blink_window #update Progress
			progress = string(floor(Int, block / Nblocks * 100))
			@js_ w (@var progress=$progress; 
					document.getElementById("simul_progress").style.width=progress+"%"; 
					document.getElementById("simul_progress").innerHTML=progress+"%"; 
					document.getElementById("simul_progress").setAttribute("aria-valuenow", progress);)
		end
	end
	#Output
	t_interp = get_sample_times(seq) 
	S_interp = LinearInterpolation(t,S,extrapolation_bc=0)(t_interp) .* get_sample_phase_compensation(seq)
	(S_interp, M0)
end

function simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams=Dict{String,Any}(), w=nothing)
	#Simulation params
	enable_gpu = get(simParams, "gpu", true) && has_cuda()
	gpu(x) = enable_gpu ? CuArray(x) : x
	Nthreads = get(simParams, "Nthreads", enable_gpu ? 1 : Hwloc.num_physical_cores())
	Δt    = get(simParams, "Δt", 1e-3)
	Δt_rf = get(simParams, "Δt_rf", 1e-4)
	t, Δt = get_uniform_times(seq, Δt; Δt_rf)
	return_Mag = get(simParams, "return_Mag", false)
	end_sim_at = get(simParams, "end_sim_at", Inf)
	if 0 < end_sim_at < dur(seq)
		idx = t .< end_sim_at
		t  =  t[idx]
		Δt = Δt[idx]
	end
	Nt = length(t)
	Nspins = prod(size(obj)...)
	Nblocks = get(simParams, "Nblocks", ceil(Int, 6506*Nt/1.15e6))
    #Simulate
	println("")
	@info "Running simulation... [GPU = $(enable_gpu), CPU = $Nthreads thread(s)]."
	if return_Mag
		_, out = @time run_sim_time_iter(obj,seq,t,Δt;Nblocks,Nthreads,gpu)
	else
		S, _ = @time run_sim_time_iter(obj,seq,t,Δt;Nblocks,Nthreads,gpu,w)
		out = S ./ Nspins #Acquired data
	end
	out
end

"""Returns magnetization of spins distributed along `z` after running the Sequence `seq`."""
function simulate_slice_profile(seq; z=range(-2e-2,2e-2,200), simParams=Dict{String,Any}("Δt_rf"=>1e-6))
	simParams["return_Mag"] = true
	sys = Scanner()
	phantom = Phantom(;x=zeros(size(z)),z)
	M = simulate(phantom, seq, sys; simParams)
	M
end