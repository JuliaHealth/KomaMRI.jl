##########################
## SIMULATION FUNCTIONS ##
##########################
gpu(x) = CUDA.has_cuda_gpu() ? CuArray(x) : x
print_gpus() = begin
	@info "$(length(devices())) CUDA capable devices."
	for d = devices()
		@info name(d)
	end
end
"""
	run_sim2D(obj,seq,t)

Simulates MRI sequence `seq` on the Phantom `obj` for time points `t`.
It calculates S(t) = ∫ ρ(x,t) exp(𝒊 ϕ(x,t)) dx.
"""
function run_sim2D_spin(obj::Phantom,seq::Sequence,t::Array{Float64,1};ϕ0::Array{Float64,1}=0.)
	𝒊 = 1im; Random.seed!(1)
	t = reshape(t,1,length(t)); Δt = t[2]-t[1]
    sz = size(obj)
    Nsz, Nt = length(sz), length(t)
	#DIFFUSION
	if !all(obj.Dλ1 .== 0) && !all(obj.Dλ2 .== 0) #No diff optimization
		ηx = sqrt.(2Δt.*obj.Dλ1).*randn(sz...,Nt)
		ηy = sqrt.(2Δt.*obj.Dλ2).*randn(sz...,Nt)
		ηxp = cumsum(ηx.*cos.(obj.Dθ).-ηy.*sin.(obj.Dθ),dims=Nsz+1)
		ηyp = cumsum(ηy.*cos.(obj.Dθ).+ηx.*sin.(obj.Dθ),dims=Nsz+1)
	else
		ηxp = 0
		ηyp = 0
	end
	#SCANNER
    Gx = get_grad(seq,1,t) |> gpu
	Gy = get_grad(seq,2,t) |> gpu 
	#SIMULATION
	ϕ0 = ϕ0 |> gpu
    xt = obj.x .+ obj.ux(obj.x,obj.y,t) .+ ηxp |> gpu
    yt = obj.y .+ obj.uy(obj.x,obj.y,t) .+ ηyp |> gpu
    if is_DAC_on(seq,t) #ACQ OPTIMIZATION
		ϕ = ϕ0 .+ (2π*γ*Δt).*cumsum(xt.*Gx.+yt.*Gy, dims=Nsz+1) 
	else
		ϕ = ϕ0 .+ (2π*γ*Δt).*sum(xt.*Gx.+yt.*Gy, dims=Nsz+1) 
	end
	#SIGNAL
	t = t		 |> gpu
	ρ = obj.ρ	 |> gpu
	Δw = obj.Δw  |> gpu
	T2 = obj.T2  |> gpu
    S = sum(ρ.*exp.(-𝒊.*(ϕ .+ Δw.*t).-t./T2 ), dims=1:Nsz)[:]
    Array(S), Array(ϕ[:,end])
end
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

"""
Implementation in multiple threads by separation the spins in N_parts.
"""
function run_sim2D_spin_parallel(obj::Phantom,seq::Sequence,t::Array{Float64,1};
	ϕ0::Array{Float64,1}=0., N_parts::Int= CUDA.has_cuda_gpu() ? 1 : Threads.nthreads())

	Nt, Ns = length(t), prod(size(obj))
	S = zeros(ComplexF64, Nt)
	
    #addprocs(N_parts)
	parts = kfoldperm(Ns, N_parts, type="ordered") 

	@threads for p ∈ parts
		aux, ϕ0[p] = run_sim2D_spin(obj[p],seq,t;ϕ0=ϕ0[p])
		S .+= aux
		aux = nothing
	end
    S, ϕ0
end

"""Divides time steps in N_parts blocks. Decreases RAM usage in long sequences."""
function run_sim2D_times_iter(obj::Phantom,seq::Sequence, t::Array{Float64,1};N_parts::Int=16)
    N, Ns = length(t), prod(size(obj))
	S = zeros(ComplexF64, N)
	ϕ0 = zeros(size(obj))

    parts = kfoldperm(N,N_parts,type="ordered")
    println("Starting simulation with Nspins=$Ns and Nt=$N")
    
	@showprogress for p ∈ parts
		S[p], ϕ0 =  run_sim2D_spin_parallel(obj, seq, t[p]; ϕ0)
	end
    S
end
