###################
# DISPLAY METHODS #
###################
# plot_seq(x::Sequence) = begin
# 	τ = dur(x)
# 	M, N = size(x.GR)
# 	T = [x.GR[1,i].T for i=1:N]
# 	t = [sum(T[1:i]) for i=1:N]
# 	t = [t[floor(Int,i/2)+1] for i=0:2*N-1]
# 	t = [0; t[1:end-1]]
# 	Δt = t[2]-t[1]
# 	idx = ["x" "y" "z"]
# 	p = [plot() for i=1:3]
# 	M0, M1, M2 = get_M0_M1_M2(sum(x))
# 	DAC = get_DAC_on(sum(x),t)
# 	for j=1:size(x[1].GR,1)
# 		p[j] = plot()
# 		T = 0
# 		M0t = [M0(t)[j] for t=t][:]/τ; 
# 		M1t = [M1(t)[j] for t=t][:]/τ^2; 
# 		M2t = [M2(t)[j] for t=t][:]/τ^3; 
# 		Gmin = minimum([g.A for g ∈ sum(x).GR[j,:]])
# 		Gmax = maximum([g.A for g ∈ sum(x).GR[j,:]])
# 		# plot!([t_k0; t_k0]*1e3,[Gmin; Gmax]*1.1e3,linewidth=0.25,color=:black)
# 		PlotlyJS.scatter(x=t[:]*1e3,y=M0t*1e3)
# 		PlotlyJS.scatter(x=t[:]*1e3,y=M1t*1e3)
# 		PlotlyJS.scatter(x=t[:]*1e3,y=M2t*1e3)
# 		for i=1:length(x)
# 			tt = t[T-Δt.<t.<T+dur(x[i])+Δt]
# 			G = get_grad(x[i],j,tt.-T)
# 						plot!(tt[:]*1e3,round.(G[:]*1e3,digits=2),legend=:none);
# 			xlabel!("Time [ms]")
# 			ylabel!("G"*idx[j]*" [mT/m]")
# 			xlims!((0,dur(sum(x))).*1e3)
# 			ylims!((Gmin,Gmax).*1.1e3)
# 			T = dur(x[i]) + T
# 		end
# 	end
# 	p
# end

plot_grads(seq::Sequence) = begin
	idx = ["Gx" "Gy" "Gz"]
	M, N = size(seq.GR)
	G = [seq.GR[j,floor(Int,i/2)+1].A for i=0:2*N-1, j=1:M]
	T = [seq.GR[1,i].T for i=1:N]
	t = [sum(T[1:i]) for i=1:N]
	t = [t[floor(Int,i/2)+1] for i=0:2*N-1]
	t = [0; t[1:end-1]]
	
	l = PlotlyJS.Layout(;title="Sequence", yaxis_title="G [mT/m]",
	    xaxis_title="t [ms]",height=300)
	p = [PlotlyJS.scatter() for j=1:M]
	for j=1:size(seq.GR,1)
		p[j] = PlotlyJS.scatter(x=t*1e3, y=G[:,j]*1e3,name=idx[j],line_shape="hv")
	end
	PlotlyJS.plot(p, l)
end
plot_grads_moments(seq::Sequence) = begin
	idx = ["Gx" "Gy" "Gz"]
	M, N = size(seq.GR)
	G = [seq.GR[j,floor(Int,i/2)+1].A for i=0:2*N-1, j=1:M]
	T = [seq.GR[1,i].T for i=1:N]
	t = [sum(T[1:i]) for i=1:N]
	t = [t[floor(Int,i/2)+1] for i=0:2*N-1]
	t = [0; t[1:end-1]]
	
	l = PlotlyJS.Layout(;title="Sequence", yaxis_title="G [mT/m]",
	    xaxis_title="t [ms]",height=300)
	p = [PlotlyJS.scatter() for j=1:M]
	for j=1:size(seq.GR,1)
		p[j] = PlotlyJS.scatter(x=t*1e3, y=G[:,j]*1e3,name=idx[j],line_shape="hv")

	end
	PlotlyJS.plot(p, l)
end
# plot_Phantom(obj::Phantom,filename::String) = begin
# 	# Phantom
# 	p1 = heatmap(obj.x*1e2,obj.y*1e2,obj.ρ,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Proton density")
# 	p2 = heatmap(obj.x*1e2,obj.y*1e2,obj.T2*1e3,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("T2 [ms]")
# 	p3 = heatmap(obj.x*1e2,obj.y*1e2,obj.Δw/2π,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("df [Hz]")
# 	P(i,j) = rotz(Dθ[i,j])[1:2,1:2]; D(i,j) = [obj.Dλ1[i,j] 0;0 obj.Dλ2[i,j]]
# 	nx = [1;0]; ny = [0;1]
# 	Dx = [nx'*P(i,j)'*D(i,j)*P(i,j)*nx for i=1:size(obj.Dλ1,1),j=1:size(obj.Dλ1,2)]
# 	Dy = [ny'*P(i,j)'*D(i,j)*P(i,j)*ny for i=1:size(obj.Dλ1,1),j=1:size(obj.Dλ1,2)]
# 	p4 = heatmap(obj.x*1e2,obj.y*1e2,Dx*1e12,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Dx [um2/s]")
# 	p5 = heatmap(obj.x*1e2,obj.y*1e2,Dy*1e12,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Dy [um2/s]")
# 	p = plot(p1,p2,p3,p4,p5,size=(1300,500),layout=@layout [a b c; d e])
# 	savefig(p,filename)
# end
# plot_sim_res(obj::Phantom,SEQ::Array{Sequence},S::Array{ComplexF64},
# 	t::Array{Float64,2},filename::String,t_k0::Float64,Ga::Float64) = begin
# 	T2m = sum(obj.T2.*obj.ρ)/sum(obj.ρ)
# 	b, n = get_bvalue(SEQ[1])
# 	println("Sequence with b-value: "*string(round(b*1e-6,digits=2))*" s/mm2") # s/mm2
# 	P(i) = rotz(obj.Dθ[i])[1:2,1:2]; D(i) = [obj.Dλ1[i] 0;0 obj.Dλ2[i]]
# 	Deq = [n'*P(i)'*D(i)*P(i)*n for i=1:prod(size(obj.Dλ1))]
# 	Dm = sum(Deq.*obj.ρ)/sum(obj.ρ) #MODIFYY

# 	p = plot_grads(SEQ,t,t_k0,Ga)
# 	p3 = plot([t_k0; t_k0]*1e3,[0; 1.1],linewidth=0.25,color=:black,label="k=0")
# 	plot!(t[:]*1e3,abs.(S),label="|S|")
# 	plot!(t[:]*1e3,exp.(-t[:]/T2m),linestyle=:dash,label="exp(-t/T2)")
# 	# exp(-b*D) <-> exp(-4*π^2*(Δ-δ/3)*q'*D*q)
# 	plot!(t[:]*1e3,exp.(-t[:]/T2m.-b*Dm).*(t[:].>=(Δ+δ)),linestyle=:dash,label="exp(-t/T2-bD)")

# 	xlabel!("Time [ms]")
# 	ylabel!("Signal [a.u]")
# 	xlims!((minimum(t),minimum(t).+dur(sum(SEQ))).*1e3)
# 	p = plot(p[1],p[2],p3,size=(800,600),layout=@layout [a ; b; c])
# 	savefig(p,filename)
# end
# plot_ksapce_trajectory(ACQ::Sequence,t::Array{Float64,2},filename::String) = begin
# 	k = get_designed_kspace(ACQ)
# 	p = plot(k[:,1],k[:,2],legend=:none,aspect_ratio=:equal)
# 	k = get_actual_kspace(ACQ,t)
# 	scatter!(k[:,1],k[:,2],legend=:none,markersize=1)
# 	xlabel!("kx [1/m]"); ylabel!("ky [1/m]")
# 	savefig(p,filename)
# end
# plot_recon(kdata::Array{ComplexF64},rec::Array{ComplexF64},
# 	Δx_pix::Float64,Δy_pix::Float64,filename::String,title::String) = begin
# 	Nx, Ny = size(rec)
# 	xr = -Δx_pix*(Nx-1)/2:Δx_pix:Δx_pix*(Nx-1)/2
# 	yr = -Δy_pix*(Ny-1)/2:Δy_pix:Δy_pix*(Ny-1)/2
# 	Wx, Wy = 1/Δx_pix, 1/Δy_pix
# 	kx = range(-Wx/2,stop=Wx/2,length=Nx)
# 	ky = range(-Wy/2,stop=Wy/2,length=Ny)
# 	p1 = heatmap(kx,ky,log.(abs.(kdata).+1),aspect_ratio=:equal,
# 				legend=:none,size=(400,400))
# 	xaxis!("kx [1/m]"); yaxis!("ky [1/m]")
# 	title!("k-space")
# 	p2 = heatmap(xr*1e2,yr*1e2,abs.(rec),aspect_ratio=:equal,
# 				legend=:none,size=(400,400))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!(title)
# 	#savefig(p1,filename*"_ksp.pdf")
# 	savefig(p2,filename*".pdf")
# end
# plot_Eq(vx,vy,S,S0,SEQ) = begin
# 	Nq, Nθ = size(SEQ)
# 	rs = range(0,stop=1,length=Nq)
# 	θs = range(0,stop=π,length=Nθ); θs = [θs[1:end-1]; π.+θs[1:end-1]];
#     Eq(vx,vy) = begin
#         Eq = [S[i,j][vx,vy] for i=1:Nq,j=1:Nθ]
#         Eq = [Eq[:,1:end-1] Eq[:,1:end-1]]
#     end
#     c = RGB{Float64}(1,1,1)
#     α = real.(abs.(S0[vx,vy]))/maximum(real.(abs.(S0)))
#     pyplot()
#     hm = heatmap(θs,rs,Eq(vx,vy)*(α<.2 ? 0 : 1),proj=:polar,aspect_ratio=:equal,
#     	legend=:none,grid=false,xticks=:none,yticks=:none,
#     	background_color_subplot=α*c)
# end
# plot_Pr(vx,vy,S,S0,SEQ) = begin
# 	Nq, Nθ = size(SEQ)
# 	rs = range(0,stop=1,length=10)*25e-6
# 	θs = range(0,stop=2π,length=10)
#     Pr(vx,vy) = begin
#         Eq = [S[i,j][vx,vy] for i=1:Nq,j=1:Nθ]
# 		b = [get_bvalue(SEQ[i,j])[1] for i=1:Nq,j=1:Nθ]
# 		n = [get_bvalue(SEQ[i,j])[2] for i=1:Nq,j=1:Nθ]
# 		B = [-b[i,j]*[n[i,j][1]^2
# 					2*n[i,j][1]*n[i,j][2]
# 					  n[i,j][2]^2] for i=1:Nq,j=1:Nθ]
# 		B = [B[i][j] for i=1:Nq*Nθ,j=1:3] #reshape
# 		y = log.(abs.(Eq))[:]
# 		# exp(-b*n'D n) <-> exp(-4*π^2*(Δ-δ/3)*q'*D*q)
# 		# min D || -b*n'D*n - log(S) ||_2 + λ ||D||_2
# 		# min D || B D - y ||_2 + λ ||D||_2 => D = (B'B+λI)^-1 B' y
# 		λ = 0
# 		Id = Matrix{Float64}(I,3,3); D = [0;0;0]
# 		for n = 0:2 #Tikhonov regularized and iterativly weighted
# 			W = n==0 ? diagm(0=>Eq[:].^2) : diagm(0=>exp.(2*B*D))
# 			D = (B'*W*B + λ*Id)^-1*B'*W*y
# 		end
# 		# Diffusion propagator
# 		D_inv = [D[1] D[2];
# 				 D[2] D[3]]^-1
# 		pr = exp.(-[([r*cos(θ) r*sin(θ)]*D_inv*[r*cos(θ);r*sin(θ)])[1] for r=rs, θ=θs])
#     end
#     c = RGB{Float64}(1,1,1)
#     α = real.(abs.(S0[vx,vy]))/maximum(real.(abs.(S0)))
#     pyplot()
#     hm = heatmap(θs,rs,Pr(vx,vy)*(α<.2 ? 0 : 1),proj=:polar,aspect_ratio=:equal,
#     	legend=:none,grid=false,xticks=:none,yticks=:none,
#     	background_color_subplot=α*c)
# end
# plot_phantom_anim(obj::Phantom, map::Array{Float64}, filename, Δxr, Δx, fps,
# 	tk=0) = begin

# 	MAX = maximum(map)
# 	gridpos = abs.(mod.(obj.x*1e2, Δxr*1e2)).*abs.(mod.(obj.y*1e2,
# 											Δxr*1e2)).≤ Δx*1e2/100*Δx*1e2/100
# 	map_grid = copy(map)
# 	map_grid[gridpos,:] = MAX.-map[gridpos,:]
# 	if fps≠0
# 		t_anim = reshape(range(0,stop=1-1/fps,length=fps),1,fps)
# 		anim = @animate for tk = t_anim
# 		    xt = obj.x.+obj.ux(obj.x,obj.y,tk);
# 			yt = obj.y.+obj.uy(obj.x,obj.y,tk)
# 		    p = scatter(xt*1e2,yt*1e2,zcolor=map_grid,marker=(:o, 2),
# 		    clim=(0,MAX),aspect_ratio=:equal,size=(380,380),
# 			label="",markerstrokewidth=0)#,colorbar_title="T2 [ms]")
# 		    xtmin = minimum(obj.x.+obj.ux(obj.x,obj.y,t_anim))
# 			xtmax = maximum(obj.x.+obj.ux(obj.x,obj.y,t_anim))
# 			ytmin = minimum(obj.y.+obj.uy(obj.x,obj.y,t_anim))
# 			ytmax = maximum(obj.y.+obj.uy(obj.x,obj.y,t_anim))
# 			xlims!(1e2*xtmin,1e2*xtmax);
# 		    ylims!(1e2*ytmin,1e2*ytmax)
# 		    xlabel!("x [cm]");ylabel!("y [cm]")
# 			title!(@sprintf("t = %1.1f s",tk))
# 		end
# 		gif(anim, filename, fps = fps)
# 	else
# 		xt = obj.x.+obj.ux(obj.x,obj.y,tk);
# 		yt = obj.y.+obj.uy(obj.x,obj.y,tk)
# 	    p = scatter(xt*1e2,yt*1e2,zcolor=map_grid,marker=(:o, 2),
# 	    clim=(0,MAX),aspect_ratio=:equal,size=(380,380),
# 		label="",markerstrokewidth=0)#,colorbar_title="T2 [ms]")
# 	    xtmin = minimum(obj.x.+obj.ux(obj.x,obj.y,tk))
# 		xtmax = maximum(obj.x.+obj.ux(obj.x,obj.y,tk))
# 		ytmin = minimum(obj.y.+obj.uy(obj.x,obj.y,tk))
# 		ytmax = maximum(obj.y.+obj.uy(obj.x,obj.y,tk))
# 		xlims!(1e2*xtmin,1e2*xtmax);
# 	    ylims!(1e2*ytmin,1e2*ytmax)
# 	    xlabel!("x [cm]");ylabel!("y [cm]")
# 		title!(@sprintf("t = %1.1f s",tk))
# 	end
# end
