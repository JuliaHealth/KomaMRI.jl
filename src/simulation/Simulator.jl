##########################
## SIMULATION FUNCTIONS ##
##########################

"""Divides a list of indices 1:N in k groups"""
function kfoldperm(N,k; type="random")
	n,r = divrem(N,k)
	b = collect(1:n:N+1)
	for i in 1:length(b)
		b[i] += i > r ? r : i-1
	end
	if type=="random"
		p = randperm(N)
	elseif type=="ordered"
		p = 1:N
	end
	return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k]]
end

#GPU related functions
gpu(x) = has_cuda() ? CuArray(x) : x
print_gpus() = begin
	println( "$(length(devices())) CUDA capable device(s)." )
	for d = devices()
		println( "  - "*name(d) )
	end
end
print_gpus_info() = begin
	@info "$(length(devices())) CUDA capable device(s)."
	for d = devices()
		@info "  - "*name(d) 
	end
end

function trapz(dx,y)
	1/2*sum( dx .* (y[:, 2:end] .+ y[:, 1:end-1]); dims=2) 
end

function cumtrapz(dx,y)
	1/2*cumsum( dx .* (y[:, 2:end] .+ y[:, 1:end-1]); dims=2)
end

function points_from_key_times(times ; dt)
	t = Float64[]
	for i = 1:length(times)-1
		if dt < times[i+1] - times[i]
			taux = collect(range(times[i],times[i+1];step=dt))
		else
			taux = [times[i], times[i+1]]
		end
		append!(t, taux)
	end
	t
end
"""
Uniform time-step calculation
"""
function get_uniform_times(seq,Δt)
	t, Δt = get_variable_times(seq; dt=Δt)
end

"""
Variable time-step calculation
"""
function get_variable_times(seq; dt=0)
	t = Float64[]
	ΔT = durs(seq) #Duration of sequence block
	T0 = cumsum([0; ΔT[1:end-1]]) #Start time of each block
	for i = 1:length(seq)
		t0 = T0[i]
		T = dur(seq[i]) #Length of block
		tf = t0+T
		if is_ADC_on(seq[i])
			t1 = t0 + seq.GR.delay[i]
			t2 = t1 + seq.GR.rise[i]
			t3 = t2 + seq.GR.T[i]
			t4 = t3 + seq.GR.fall[i]
			taux = points_from_key_times([t0,t1,t2,t3,t4,tf]; dt)
		else
			t1 = t0 + seq.GR.delay[i]
			t2 = t1 + seq.GR.rise[i]
			t3 = t2 + seq.GR.T[i]
			t4 = t3 + seq.GR.fall[i]
			taux = points_from_key_times([t0,t1,t2,t3,t4,tf]; dt)
		end
		append!(t, taux)
	end
	t = unique(x -> round(x*1e8), t)
	Δt = t[2:end] .- t[1:end-1]
	t = t[1:end-1] .+ 1e-8
	t, Δt
end

"""
Implementation in multiple threads. Separating the spins in N_parts.
"""
function run_spin_precession_parallel(obj::Phantom,seq::Sequence, t::Array{Float64,1}, Δt::Array{Float64,1};
	M0::Array{Mag,1}, Nthreads::Int=Threads.nthreads())

	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	#Put times as row vector
	t = reshape(t,1,Nt)
	Δt = reshape(Δt,1,NΔt)

	S = zeros(ComplexF64, Nt)
	
	parts = kfoldperm(Ns, Nthreads, type="ordered") 

	@threads for p ∈ parts
		aux, M0[p] = run_spin_precession(obj[p],seq,t,Δt; M0=M0[p])
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
	M0::Array{Mag,1})

	𝒊 = 1im; Random.seed!(1)
	T = sum(Δt) #Total length, used for signal relaxation
	# t = [t t[end]+Δt[end]] #For trapezoidal integration <---------------------------
    sz = size(obj)
	Nsz = length(sz)

	#DIFFUSION, disabled while I think how to do it efficiently
	# if any(obj.Dλ1 .+ obj.Dλ2 .!= 0)  #No diff optimization
	# 	#TODO: I need to add diff displacement η story between blocks (like ϕ0)
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
	tp = t .- t[1]	|> gpu # t' = t - t0
	t = t			|> gpu
	Δt = Δt			|> gpu
    xt = x0 .+ obj.ux(x0,y0,z0,t) .+ ηxp |> gpu
	yt = y0 .+ obj.uy(x0,y0,z0,t) .+ ηyp |> gpu
	zt = z0 .+ obj.uz(x0,y0,z0,t) .+ ηzp |> gpu
	#ACQ OPTIMIZATION
    if is_ADC_on(seq, Array(t)) 
		ϕ = ϕ0 .- (2π*γ).*cumsum(Δt.*(xt.*Gx .+ yt.*Gy .+ zt.*Gz), dims=2)
	else
		ϕ = ϕ0 .- (2π*γ).*sum(Δt.*(xt.*Gx .+ yt.*Gy .+ zt.*Gz), dims=2) 
	end
	#Mxy preccesion and relaxation
	Δw = obj.Δw  |> gpu
	T2 = obj.T2  |> gpu
	Mxy = Mxy .* exp.(𝒊.*(ϕ .- Δw.*tp) .- tp./T2 )
	#ACQUIRED SIGNAL
	S = sum(Mxy, dims=1:Nsz)[:] #<--- TODO: add coil sensitivities
	#Mz relaxation
	T1 = obj.T1		|> gpu
	Mz0 = obj.ρ		|> gpu
	Mz = M0.z		|> gpu
	Mz = Mz .* exp.(-T./T1) .+ Mz0 .* ( 1 .- exp.(-T./T1) )
	# END
	Mxy = Array(Mxy)
	Mz = Array(Mz)
	M0 = Mag.(Mxy[:,end], Mz) #Saving the last magnetization
    Array(S), M0
end

run_spin_excitation_parallel(obj, seq, t::Array{Float64,1}, Δt::Array{Float64,1}; 
	M0::Array{Mag,1}, Nthreads::Int=Threads.nthreads()) = begin
	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	#Put times as row vector
	t = reshape(t,1,Nt)
	Δt = reshape(Δt,1,NΔt)
	
	parts = kfoldperm(Ns, Nthreads, type="ordered") 

	@threads for p ∈ parts
		M0[p] = run_spin_excitation(obj[p],seq,t,Δt; M0=M0[p])
	end
    M0
end

run_spin_excitation(obj, seq, t::Array{Float64,2}, Δt::Array{Float64,2}; 
	M0::Array{Mag,1}) = begin
	#SCANNER
	B1 = 		get_rfs(seq,t)[1]
    Gx, Gy, Gz = 	get_grads(seq,t)
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
	ΔB0 = obj.Δw./(2π*γ)				|> gpu
	Bz = (Gx.*xt .+ Gy.*yt .+ Gz.*zt) .+ ΔB0	#<-- This line is very slow, FIX!!
	B = sqrt.(abs.(B1).^2. .+ abs.(Bz).^2.)		
	φ = -2π*γ * (B .* Δt) # angle of rotation 
	B[B.==0] .= 1e-17; # removes problems when dividing by φ
	Qt = Q.(Array(φ), Array(B1./B), Array(Bz./B))
	Qf = prod( Qt , dims=2 )[:] # equivalent rotation
	#TODO: Relaxation effects
	M0 =  Qf .* M0
end

"""Divides time steps in N_parts blocks. Decreases RAM usage in long sequences."""
function run_sim_time_iter(obj::Phantom,seq::Sequence, t::Array{Float64,1}, Δt; 
								Nblocks::Int=16, Nthreads::Int=Threads.nthreads())
	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	if NΔt ==1 Δt = Δt*ones(size(t)) end
	#Put times as row vector
	S = zeros(ComplexF64, Nt)
	if is_RF_on(seq)
		M0 = Mag(obj,:z)
	else
		M0 = Mag(obj,:x)
	end
    parts = kfoldperm(Nt,Nblocks,type="ordered")
	println("Starting simulation with Nspins=$Ns and Nt=$Nt")
	#TODO: transform suceptibility χ to Δω, for each time-block with FMM-like technique O(nlogn).
	@showprogress for p ∈ parts
		if is_RF_on(seq, t[p])
			M0  = run_spin_excitation_parallel(obj, seq, t[p], Δt[p]; M0, Nthreads)
		else
			S[p], M0 = run_spin_precession_parallel(obj, seq, t[p], Δt[p]; M0, Nthreads)
		end
	end

	#TODO: output raw data in ISMRMD format
	t_interp = get_sample_times(seq)
	S_interp = LinearInterpolation(t.+Δt,S)(t_interp)
	(S_interp, t_interp)
end

function simulate(phantom::Phantom, seq::Sequence, simParams::Dict)
	#Simulation params
	Nthreads = get(simParams, "Nthreads", Threads.nthreads())
	step = get(simParams, "step", "uniform")
	if step == "uniform"
		Δt = get(simParams, "Δt", 4e-6) #<- simulate param
		t, Δt = get_uniform_times(seq,Δt)
		Nspins, Nt = prod(size(phantom)), length(t)
		Nblocks = get(simParams, "Nblocks", floor(Int, Nspins*Nt/2.7e6))
	elseif step == "variable"
		t, Δt = get_variable_times(seq)
		Nblocks = floor(Int64, length(t) / 1)
	end
	println("Dividing simulation in Nblocks=$Nblocks")
    #Simulate
    S, t_interp = @time run_sim_time_iter(phantom,seq,t,Δt;Nblocks,Nthreads)
    Nspins = prod(size(phantom))
	signal = S ./ Nspins #Acquired data
	#Recon params parsing
	recParams = Dict()
    recParams["Nx"] =  get(seq.DEF, "Nx",  101)
	recParams["Ny"] =  get(seq.DEF, "Ny",  recParams["Nx"])
	recParams["EPI"] = get(seq.DEF, "EPI", false)
	(signal, t_interp, recParams)
end