##########################
## SIMULATION FUNCTIONS ##
##########################
# @everywhere begin
"""
	run_sim2D(obj,seq,t)

Simulates MRI sequence `seq` on the Phantom `obj` for time points `t`.
It calculates S(t) = ∫ ρ(x,t) exp(𝒊 ϕ(x,t)) dx.
"""
run_sim2D_spin(obj::Phantom,Seq::Sequence,t::Array{Float64,2};ϕ0::Array{Float64,1}=0.) = begin
	𝒊 = 1im; Random.seed!(1) # Setting random seed for comparable results
	sz = size(obj);	Nsz = length(sz)
	Δt = t[2]-t[1]; Nt = length(t)

	#DIFF MODEL
	# Initial position + Diffusion + Displacement field
	ηx = sqrt.(2Δt.*obj.Dλ1).*randn(sz...,Nt)
	ηy = sqrt.(2Δt.*obj.Dλ2).*randn(sz...,Nt)
	ηxp = cumsum(ηx.*cos.(obj.Dθ).-ηy.*sin.(obj.Dθ),dims=Nsz+1)
	ηyp = cumsum(ηy.*cos.(obj.Dθ).+ηx.*sin.(obj.Dθ),dims=Nsz+1)

	#SCANNER
	Gx = get_grad(Seq,1,t) #<-scanner
	Gy = get_grad(Seq,2,t) #<-scanner

	#SIMULATION
	xt = obj.x .+ obj.ux(obj.x,obj.y,t) .+ ηxp
	yt = obj.y .+ obj.uy(obj.x,obj.y,t) .+ ηyp
	ϕ = ϕ0 .+ (2π*γ*Δt).*cumsum(xt.*Gx.+yt.*Gy, dims=Nsz+1) #SLOW!
	S = sum(obj.ρ.*exp.(-𝒊.*(ϕ .+ obj.Δw.*t).-t.*obj.T2.^-1 ), dims=1:Nsz)[:] #MRI-signal with T2

	#Signal, final-phase
	S, ϕ[:,end]
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
Work in progress. Implementation in multiple threads by separation the spins in N_parts.
"""
run_sim2D_spin_parallel(obj::Phantom,Seq::Sequence,
	t::Array{Float64,2},N_parts::Int=4) = begin
	S = zeros(ComplexF64,length(t))
	N = length(obj.ρ)
	
	addprocs(N_parts)
	parts = kfoldperm(N, N_parts)
	@everywhere sub_part(obj::Phantom,p::Array) = begin
		Phantom(obj.name,obj.x[p],obj.y[p],obj.ρ[p],
					obj.T2[p],obj.Δw[p],obj.Dλ1[p],obj.Dλ2[p],
					obj.Dθ[p],obj.ux,obj.uy)
	end
	S = @distributed (+) for p ∈ parts
		run_sim2D_spin(sub_part(obj,p),Seq,t)
	end
	S
end

"""Divides time steps in N_parts blocks, to decreease RAM usage in long sequences."""
run_sim2D_times_iter(obj::Phantom,Seq::Sequence,
	t::Array{Float64,2},N_parts::Int=8) = begin
	N, Ns = length(t), prod(size(obj))
	S = zeros(ComplexF64, N)
	parts = kfoldperm(N,N_parts,type="ordered")
	
	println("Starting simulation with Nspins=$Ns and Nt=$N")
	ϕ0 = zeros(size(obj)) #phase history
	@showprogress for p ∈ parts
		S[p], ϕ0 = run_sim2D_spin(obj,Seq, reshape(t[p],1,length(p)), ϕ0=ϕ0)
	end
	S
end
