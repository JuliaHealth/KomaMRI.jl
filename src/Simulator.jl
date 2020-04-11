##########################
## SIMULATION FUNCTIONS ##
##########################
# @everywhere begin
"""
	run_sim2D(obj,seq,t)

Simulates MRI sequence `seq` on the Phantom `obj` for time points `t`.
It calculates S(t) = ∫ ρ(x,t) exp(𝒊 ϕ(x,t)) dx.
"""
run_sim2D_spin(obj::Phantom,Seq::Sequence,t::Array{Float64,2}) = begin
	𝒊 = 1im; Random.seed!(1) # Setting random seed for comparable results
	sz = size(obj);	Nsz = length(sz)
	Δt = t[2]-t[1]; Nt = length(t)
	println("####################################################################")
	println("Starting simulation for Ns="*string(prod(sz))*" spins and Nt="*string(Nt)*" time steps...")
	# Initial position + Diffusion + Displacement field
	ηx = @time sqrt.(2Δt.*obj.Dλ1).*randn(sz...,Nt)
	ηy = @time sqrt.(2Δt.*obj.Dλ2).*randn(sz...,Nt)
	Gx = get_grad(Seq,1,t)
	Gy = get_grad(Seq,2,t)
	# SLOW ->
	xt = @time obj.x.+cumsum(ηx.*cos.(obj.Dθ).-ηy.*sin.(obj.Dθ),dims=Nsz+1).+obj.ux(obj.x,obj.y,t)
	yt = @time obj.y.+cumsum(ηy.*cos.(obj.Dθ).+ηx.*sin.(obj.Dθ),dims=Nsz+1).+obj.uy(obj.x,obj.y,t)
	ϕ = @time (2π*γ*Δt).*cumsum(xt.*Gx.+yt.*Gy, dims=Nsz+1) #SLOW!
	S = @time sum(obj.ρ.*exp.(-𝒊.*(ϕ.+obj.Δw.*t).-t.*obj.T2.^-1 ), dims=1:Nsz)[:] #MRI-signal with T2
	println("Simulation completed!")
	println("####################################################################")
	print("Total simulation time: ")
	S
end
# end
function kfoldperm(N,k)
	n,r = divrem(N,k)
	b = collect(1:n:N+1)
	for i in 1:length(b)
		b[i] += i > r ? r : i-1
	end
	p = randperm(N)
	return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k]]
end
run_sim2D_spin_parallel(obj::Phantom,Seq::Sequence,
	t::Array{Float64,2},N_parts::Int=4) = begin
	S = zeros(ComplexF64,length(t))
	N = length(obj.ρ)
	addprocs(N_parts)
	parts = kfoldperm(N,N_parts)
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
