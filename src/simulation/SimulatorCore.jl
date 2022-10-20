##########################
## SIMULATION FUNCTIONS ##
##########################

#GPU related functions
"""
    print_gpus()

Simple function to print the gpus CUDA devices available in the host.
"""
print_gpus() = begin
	println( "$(length(devices())) CUDA capable device(s)." )
	for (i,d) = enumerate(devices())
		u = i == 1 ? "*" : " "
		println( "  ($(i-1)$u) $(name(d))")
	end
end

"""
    S, M0 = run_spin_precession_parallel(obj, seq, t, Δt; M0, Nthreads, gpu)

Implementation in multiple threads for the simulation in free precession,
separating the spins of the phantom `obj` in `Nthreads`.

!!! note
    This function does not use the complete time vector of the total simulation, it uses
    instead a piece of the time vector given by a variable called `Nblocks` (which is
    outside of this function, refer to [`run_sim_time_iter`](@ref)) to reduce the RAM memory
    utilization.

# Arguments
- `obj`: (`::Phantom`) the phantom struct
- `seq`: (`::Sequence`) the sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) the non-uniform time vector (actually it's a part of the
    complete simulation time vector)
- `Δt`: (`::Vector{Float64}`, `[s]`) the delta time of `t` (actually it's a part of the
    complete simulation time vector)

# Keywords
- `M0`: (`::Vector{Mag}`) the initial state of the Mag vector
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) the number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) the function that represents the gpu of the host

# Returns
- `S`: (`Vector{ComplexF64}`) the raw signal over time
- `M0`: (`::Vector{Mag}`) the final state of the Mag vector (or the initial state for the
    next simulation step (the next step can be another precession step or an excitation
    step))
"""
function run_spin_precession_parallel(obj::Phantom, seq::Sequence, t::Array{Float64,1}, Δt::Array{Float64,1};
	M0::Array{Mag,1}, Nthreads::Int=Hwloc.num_physical_cores(), gpu::Function)

	Nt, NΔt, Ns = length(t), length(Δt), prod(size(obj))
	#Put times as row vector
	t = reshape(t,1,Nt)
	Δt = reshape(Δt,1,NΔt)

	S = zeros(ComplexF64, Nt)

	parts = kfoldperm(Ns, Nthreads, type="ordered")

	S = ThreadsX.mapreduce(+, parts) do p #Thread-safe summation
		S_p, M0[p] = run_spin_precession(obj[p],seq,t,Δt; M0=M0[p], gpu)
		S_p
	end

    S, M0
end

"""
    S, M0 = run_spin_precession(obj, seq, t, Δt; M0, gpu)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∫ ρ(x,t) exp(- t/T2(x,t) ) exp(- 𝒊 ϕ(x,t)) dx. It performs the simulation in free
precession.

!!! note
    This function is used to simulate a part of the simulation over time given by the
    variable `Nblocks` (which is outside of this function, refer to
    [`run_sim_time_iter`](@ref)) to reduce the RAM memory utilization. It is also used to
    simulate a part of the spins in a phantom defined by the variable `Nthreads` (which is
    outside of this function too, refer to [`run_spin_precession_parallel`](@ref)) to take
    advantage of CPU parallel processing.

# Arguments
- `obj`: (`::Phantom`) the phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) the sequence struct
- `t`: (`1-row ::Matrix{Float64}`, `[s]`) the non-uniform time vector (actually it's a part
    of the complete simulation time vector)
- `Δt`: (`1-row ::Matrix{Float64}`, `[s]`) the delta time of `t` (actually it's a part of
    the complete simulation time vector)

# Keywords
- `M0`: (`::Vector{Mag}`) the initial state of the Mag vector (actually, it's a part of the
    complete Mag vector)
- `gpu`: (`::Function`) the function that represents the gpu of the host

# Returns
- `S`: (`Vector{ComplexF64}`) the raw signal over time
- `M0`: (`::Vector{Mag}`) the final state of the Mag vector (actually, it's a part of the
    complete Mag vector) (it's not the initial state for the next simulation, since it's
    necessary to add the magnetization of all the parts of the phantom (i.e. sum up all the
    spin magnetizations first), refer to [`run_spin_precession_parallel`](@ref))
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

"""
    M0 = run_spin_excitation_parallel(obj, seq, t, Δt; M0, Nthreads, gpu)

It gives rise to a rotation of M0 with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis. It uses
different number threads to excecute the process.

!!! note
    This function does not use the complete time vector of the total simulation, it uses
    instead a piece of the time vector given by a variable called `Nblocks` (which is
    outside of this function, refer to [`run_sim_time_iter`](@ref)) to reduce the RAM memory
    utilization.

# Arguments
- `obj`: (`::Phantom`) the phantom struct
- `seq`: (`::Sequence`) the sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) the non-uniform time vector (actually it's a part of the
    complete simulation time vector)
- `Δt`: (`::Vector{Float64}`, `[s]`) the delta time of `t` (actually it's a part of the
    complete simulation time vector)

# Keywords
- `M0`: (`::Vector{Mag}`) the initial state of the Mag vector
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) the number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) the function that represents the gpu of the host

# Returns
- `M0`: (`::Vector{Mag}`) the final state of the Mag vector after a rotation (or the initial
    state for the next precession simulation step)
"""
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

"""
    M0 = run_spin_excitation(obj, seq, t, Δt; M0, gpu)

It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis.

!!! note
    This function is used to simulate a part of the simulation over time given by the
    variable `Nblocks` (which is outside of this function, refer to
    [`run_sim_time_iter`](@ref)) to reduce the RAM memory utilization. It is also used to
    simulate a part of the spins in a phantom defined by the variable `Nthreads` (which is
    outside of this function too, refer to [`run_spin_excitation_parallel`](@ref)) to take
    advantage of CPU parallel processing.

# Arguments
- `obj`: (`::Phantom`) the phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) the sequence struct
- `t`: (`1-row ::Matrix{Float64}`, `[s]`) the non-uniform time vector (actually it's a part
    of the complete simulation time vector)
- `Δt`: (`1-row ::Matrix{Float64}`, `[s]`) the delta time of `t` (actually it's a part of
    the complete simulation time vector)

# Keywords
- `M0`: (`::Vector{Mag}`) the initial state of the Mag vector (actually, it's a part of the
    complete Mag vector)
- `gpu`: (`::Function`) the function that represents the gpu of the host

# Returns
- `M0`: (`::Vector{Mag}`) the final state of the Mag vector after a rotation (actually, it's
    a part of the complete Mag vector and it's a part of the initial state for the next
    precession simulation step)
"""
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

"""
    S_interp, M0 = run_sim_time_iter(obj, seq, t, Δt; Nblocks, Nthreads, gpu, w)

Performs the simulation over the total time vector `t` by dividing the time into `Nblocks`
parts to reduce RAM usage and spliting the spins of the phantom `obj` into `Nthreads` to
take advantage of CPU parallel processing.

# Arguments
- `obj`: (`::Phantom`) the phantom struct
- `seq`: (`::Sequence`) the sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) the non-uniform time vector
- `Δt`: (`::Vector{Float64}`, `[s]`) the delta time of `t`

# Keywords
- `Nblocks`: (`::Int`, `=16`) the number of groups for spliting the simulation over time
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) the number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) the function that represents the gpu of the host
- `w`: (`::Any`, `=nothing`) the flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `S_interp`: (`::Vector{ComplexF64}`) the interpolated raw signal
- `M0`: (`::Vector{Mag}`) the final state of the Mag vector
"""
function run_sim_time_iter(obj::Phantom, seq::Sequence, t::Array{Float64,1}, Δt;
								Nblocks::Int=16, Nthreads::Int=Hwloc.num_physical_cores(), gpu::Function, w=nothing)
	Nt, Ns = length(t), prod(size(obj))
	blink_window = w !== nothing
	S = zeros(ComplexF64, Nt) #Only one coil for now, TODO: change for more coils to Nt x Coils
	M0 = Mag(obj,:z)
	breaks = get_breaks_in_RF_key_points(seq,t)
    parts = kfoldperm(Nt,Nblocks;type="ordered",breaks)
	Nblocks = length(parts)
	# To visually check the simulation blocks
	t_sim_parts = [t[p[1]] for p in parts]
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
	S_interp = LinearInterpolation(t.+Δt,S,extrapolation_bc=0)(t_interp) .* get_sample_phase_compensation(seq)
	(S_interp, M0, t_sim_parts)
end

"""
    out = simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams, w)

Returns the raw signal or the last state of the magnetization according to the value
of the `"return_type"` key of the `simParams` dictionary.

# Arguments
- `obj`: (`::Phantom`) the phantom struct
- `seq`: (`::Sequence`) the sequence struct
- `sys`: (`::Scanner`) the scanner struct

# Keywords
- `simParams`: (`::Dict{String,Any}`, `=Dict{String,Any}()`) the dictionary with simulation
    parameters
- `w`: (`::Any`, `=nothing`) the flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `out`: (`::Vector{ComplexF64}` or `::Vector{Mag}` or `RawAcquisitionData`) depending if
    "return_type" is "mat" or "mag" or "raw" (default) respectively.

# Examples

Preparation (define scanner and sequence):
```julia-repl
julia> sys = Scanner();

julia> FOV, N = 23e-2, 101;

julia> durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse

julia> ex = PulseDesigner.RF_hard(sys.B1, durRF, sys)
Sequence[ τ = 0.587 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> epi = PulseDesigner.EPI(FOV, N, sys)
Sequence[ τ = 62.259 ms | blocks: 203 | ADC: 101 | GR: 205 | RF: 0 | DEF: 4 ]

julia> seq = ex + epi
Sequence[ τ = 62.846 ms | blocks: 204 | ADC: 101 | GR: 205 | RF: 1 | DEF: 4 ]

julia> plot_seq(seq)

julia> plot_kspace(seq)
```

Simulate:
```julia-repl
julia> obj = brain_phantom2D()

julia> ismrmrd = simulate(obj, seq, sys);

julia> plot_signal(ismrmrd)
```

Reconstruct:
```julia-repl
julia> Nx, Ny = ismrmrd.params["reconSize"][1:2];

julia> params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true);

julia> acq = AcquisitionData(ismrmrd);

julia> recon = reconstruction(acq, params);

julia> image = reshape(recon.data, Nx, Ny, :)
102×102×1 Array{ComplexF64, 3}:
[:, :, 1] =
 0.0+0.0im  0.0+0.0im  …  0.0+0.0im
 0.0+0.0im  0.0+0.0im     0.0+0.0im
    ⋮           ⋮       ⋱      ⋮
 0.0+0.0im  0.0+0.0im  …  0.0+0.0im

julia> slice_abs = abs.(image[:, :, 1])
102×102 Matrix{Float64}:
 0.0  0.0  …  0.0
 0.0  0.0     0.0
  ⋮        ⋱   ⋮
 0.0  0.0  …  0.0

julia> plot_image(slice_abs)
```
```
"""
function simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams=Dict{String,Any}(), w=nothing)
	#Simulation params
	enable_gpu = get(simParams, "gpu", true) && has_cuda()
	gpu(x) = enable_gpu ? CuArray(x) : x
	Nthreads = get(simParams, "Nthreads", enable_gpu ? 1 : Hwloc.num_physical_cores())
	Δt    = get(simParams, "Δt", 1e-3)
	Δt_rf = get(simParams, "Δt_rf", 1e-4)
	t, Δt = get_uniform_times(seq, Δt; Δt_rf)
	return_type = get(simParams, "return_type", "raw")
	end_sim_at =  get(simParams, "end_sim_at", Inf)
	if 0 < end_sim_at < dur(seq)
		idx = t .< end_sim_at
		t  =  t[idx]
		Δt = Δt[idx]
	end
	Nt = length(t)
	Nspins = prod(size(obj)...)
	Nblocks = get(simParams, "Nblocks", ceil(Int, 6506*Nt/1.15e6))
    #Simulate
	@info "Running simulation... [GPU = $(enable_gpu), CPU = $Nthreads thread(s)]."
	@time begin
		timed_tuple = @timed run_sim_time_iter(obj,seq,t,Δt;Nblocks,Nthreads,gpu,w)
	end
	S, M, t_sim_parts = timed_tuple.value #unpacking
	out = S ./ Nspins #Acquired data
	if return_type == "mag"
		out = M
	elseif return_type == "mat"
		out = S
	elseif return_type == "raw"
		simParams_raw = copy(simParams)
		simParams_raw["gpu"] = enable_gpu
		simParams_raw["Nthreads"] = Nthreads
		simParams_raw["t_sim_parts"] = t_sim_parts
		simParams_raw["Nblocks"] = Nblocks
		simParams_raw["sim_time"] = timed_tuple.time
		out = signal_to_raw_data([S;;], seq; phantom=obj, sys=sys, simParams=simParams_raw)
	end
	out
end

"""
    M = simulate_slice_profile(seq; z, simParams)

Returns magnetization of spins distributed along `z` after running the Sequence `seq`.

!!! note
    This function is not being used in this KomaMRI version.

# Arguments
- `seq`: (`::Sequence`) the sequence struct

# Keywords
- `z`: (`=range(-2e-2,2e-2,200)`) a range for the z axe
- `simParams`: (`::Dict{String, Any}`, `=Dict{String,Any}("Δt_rf"=>1e-6)`) a dictionary with
    simulation parameters

# Returns
- `M`: (`::Vector{Mag}`) the final state of the Mag vector
"""
function simulate_slice_profile(seq; z=range(-2e-2,2e-2,200), simParams=Dict{String,Any}("Δt_rf"=>1e-6))
	simParams["return_type"] = "raw"
	sys = Scanner()
	phantom = Phantom(;x=zeros(size(z)),z)
	M = simulate(phantom, seq, sys; simParams)
	M
end
