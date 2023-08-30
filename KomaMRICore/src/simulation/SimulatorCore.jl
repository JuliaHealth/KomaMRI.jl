abstract type SimulationMethod end #get all available types by using subtypes(KomaMRI.SimulationMethod)
abstract type SpinStateRepresentation{T<:Real} end #get all available types by using subtypes(KomaMRI.SpinStateRepresentation)

#Defined methods:
include("Bloch/BlochSimulationMethod.jl")       #Defines Bloch simulation method
include("Bloch/BlochDictSimulationMethod.jl")   #Defines BlochDict simulation method
include("Bloch/BlochMovSimulationMethod.jl")    #Defines BlochMov simulation method

"""
    sig, Xt = run_spin_precession_parallel(obj, seq, M; Nthreads)

Implementation in multiple threads for the simulation in free precession,
separating the spins of the phantom `obj` in `Nthreads`.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `sig`: (`Vector{ComplexF64}`) raw signal over time
- `Xt`: (`::Vector{Mag}`) final state of the Mag vector (or the initial state for the
    next simulation step (the next step can be another precession step or an excitation
    step))
"""
function run_spin_precession_parallel!(obj::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod,
    Ux,Uy,Uz;
    Nthreads = Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads, type="ordered")
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        ux = Ux !== nothing ? @view(Ux[p,:]) : nothing
        uy = Uy !== nothing ? @view(Uy[p,:]) : nothing
        uz = Uz !== nothing ? @view(Uz[p,:]) : nothing

        run_spin_precession!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method, 
                             ux,uy,uz)
    end

    return nothing
end

"""
    M0 = run_spin_excitation_parallel(obj, seq, Xt; Nthreads)

It gives rise to a rotation of M0 with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis. It uses
different number threads to excecute the process.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (or the initial
    state for the next precession simulation step)
"""
function run_spin_excitation_parallel!(obj::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod,
    Ux,Uy,Uz;
    Nthreads=Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads; type="ordered")
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        ux = Ux !== nothing ? @view(Ux[p,:]) : nothing
        uy = Uy !== nothing ? @view(Uy[p,:]) : nothing
        uz = Uz !== nothing ? @view(Uz[p,:]) : nothing

        run_spin_excitation!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method, 
                             ux,uy,uz)
    end

    return nothing
end

"""
    S_interp, M0 = run_sim_time_iter(obj, seq, t, Δt; Nblocks, Nthreads, gpu, w)

Performs the simulation over the total time vector `t` by dividing the time into `Nblocks`
parts to reduce RAM usage and spliting the spins of the phantom `obj` into `Nthreads` to
take advantage of CPU parallel processing.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) non-uniform time vector
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time of `t`

# Keywords
- `Nblocks`: (`::Int`, `=16`) number of groups for spliting the simulation over time
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) function that represents the gpu of the host
- `w`: (`::Any`, `=nothing`) flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `S_interp`: (`::Vector{ComplexF64}`) interpolated raw signal
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_sim_time_iter!(obj::Phantom, seq::DiscreteSequence, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod,itp;
    Nblocks=1, Nthreads=Threads.nthreads(), parts=[1:length(seq)], w=nothing) where {T<:Real}
    # Simulation
    rfs = 0
    samples = 1
    progress_bar = Progress(Nblocks)

    Ux,Uy,Uz = get_displacements_2(obj,seq.t,itp)

    for (block, p) = enumerate(parts) 
        seq_block = @view seq[p]
        ux = (Ux !== nothing) ? KomaMRICore.CUDA.hcat(@view(Ux[:,p]),@view(Ux[:,p[end]])) : nothing # We need to duplicate the last column of Ux, Uy and Uz
        uy = (Uy !== nothing) ? KomaMRICore.CUDA.hcat(@view(Uy[:,p]),@view(Uy[:,p[end]])) : nothing # so that dimensions match
        uz = (Uz !== nothing) ? KomaMRICore.CUDA.hcat(@view(Uz[:,p]),@view(Uz[:,p[end]])) : nothing

        # Params
        excitation_bool = is_RF_on(seq_block) && is_ADC_off(seq_block) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
        Nadc = sum(seq_block.ADC)
        acq_samples = samples:samples+Nadc-1
        dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times
        # Simulation wrappers
        if excitation_bool
            run_spin_excitation_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method,
                                          ux,uy,uz;
                                          Nthreads)
            rfs += 1
        else
            run_spin_precession_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method,
                                          ux,uy,uz; 
                                          Nthreads)
        end
        samples += Nadc
        #Update progress
        next!(progress_bar, showvalues=[(:simulated_blocks, block), (:rf_blocks, rfs), (:acq_samples, samples-1)])
        update_blink_window_progress!(w, block, Nblocks)
    end
    return nothing
end

"""Updates KomaUI's simulation progress bar."""
function update_blink_window_progress!(w::Nothing, block, Nblocks)
    return nothing
end

"""
    out = simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams, w)

Returns the raw signal or the last state of the magnetization according to the value
of the `"return_type"` key of the `simParams` dictionary.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `simParams`: (`::Dict{String,Any}`, `=Dict{String,Any}()`) the dictionary with simulation
    parameters
- `w`: (`::Any`, `=nothing`) the flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `out`: (`::Vector{ComplexF64}` or `::S <: SpinStateRepresentation` or `RawAcquisitionData`) depending if
    "return_type" is "mat" or "state" or "raw" (default) respectively.

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq");

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> raw = simulate(obj, seq, sys)

julia> plot_signal(raw)
```
"""

function simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams=Dict{String,Any}(), w=nothing)
    #Simulation parameter parsing, and setting defaults
    enable_gpu  = get(simParams, "gpu", true); if enable_gpu check_use_cuda(); enable_gpu &= use_cuda[] end
    gpu_device  = get(simParams, "gpu_device", 0)
    Nthreads    = get(simParams, "Nthreads", enable_gpu ? 1 : Threads.nthreads())
    Nblocks     = get(simParams, "Nblocks", 20)
    Δt          = get(simParams, "Δt", 1e-3)
    Δt_rf       = get(simParams, "Δt_rf", 5e-5)
    sim_method  = get(simParams, "sim_method", Bloch())
    precision   = get(simParams, "precision", "f32")
    return_type = get(simParams, "return_type", "raw")
    # Simulation init
    t, Δt = get_uniform_times(seq, Δt; Δt_rf)
    t = [t; t[end]+Δt[end]]
    breaks = get_breaks_in_RF_key_points(seq, t)
    parts = kfoldperm(length(Δt), Nblocks; type="ordered", breaks)
    Nblocks = length(parts)
    t_sim_parts = [t[p[1]] for p ∈ parts]
    # Sequence init
    B1, Δf     = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    ADCflag = [any(tt .== tadc) for tt in t[2:end]] #Displaced 1 dt, sig[i]=S(ti+dt)
    seqd = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    # Spins' state init (Magnetization, EPG, etc.), could include modifications to obj (e.g. T2*)
    Xt, obj = initialize_spins_state(obj, sim_method)
    # Signal init
    Ndims = sim_output_dim(obj, seq, sys, sim_method)
    sig = zeros(ComplexF64, Ndims..., Nthreads)

    # NEW <------- INTERPOLATION MOTION FUNCTIONS -------------
    Ns = length(obj.x)
    limits = get_pieces_limits_cpu(obj)

    Δ = zeros(Ns,length(limits),3)
    Δ[:,:,1] = hcat(zeros(Ns,1),obj.Δx,zeros(Ns,1))
    Δ[:,:,2] = hcat(zeros(Ns,1),obj.Δy,zeros(Ns,1))
    Δ[:,:,3] = hcat(zeros(Ns,1),obj.Δz,zeros(Ns,1))

    itpx = reshape(sum(abs.(Δ[:,:,1]);dims=2),(Ns,)) != zeros(Ns) ? [interpolate((limits,), Δ[i,:,1], Gridded(Linear())) for i in 1:Ns] : nothing
    itpy = reshape(sum(abs.(Δ[:,:,2]);dims=2),(Ns,)) != zeros(Ns) ? [interpolate((limits,), Δ[i,:,2], Gridded(Linear())) for i in 1:Ns] : nothing
    itpz = reshape(sum(abs.(Δ[:,:,3]);dims=2),(Ns,)) != zeros(Ns) ? [interpolate((limits,), Δ[i,:,3], Gridded(Linear())) for i in 1:Ns] : nothing

    itp = [itpx, itpy, itpz]
    
    # --------------------------------------------------------- 
    # Precision
    if precision == "f32" #Default
        obj  = obj  |> f32 #Phantom
        seqd = seqd |> f32 #DiscreteSequence
        Xt   = Xt   |> f32 #SpinStateRepresentation
        sig  = sig  |> f32 #Signal
        # NEW <----------------------
        itp  = itp  |> f32 #Motion
        # ---------------------------
    elseif precision == "f64"
        obj  = obj  |> f64 #Phantom
        seqd = seqd |> f64 #DiscreteSequence
        Xt   = Xt   |> f64 #SpinStateRepresentation
        sig  = sig  |> f64 #Signal
        # NEW <----------------------
        itp  = itp  |> f64 #Motion
        # ---------------------------
    end
    # Objects to GPU
    if enable_gpu #Default
        device!(gpu_device)
        gpu_name = name.(devices())[gpu_device+1]
        obj  = obj  |> gpu #Phantom
        seqd = seqd |> gpu #DiscreteSequence
        Xt   = Xt   |> gpu #SpinStateRepresentation
        sig  = sig  |> gpu #Signal
        # NEW <----------------------
        itp  = itp  |> gpu; #Motion
        # ---------------------------
    end

    # Simulation
    @info "Running simulation in the $(enable_gpu ? "GPU ($gpu_name)" : "CPU with $Nthreads thread(s)")" koma_version=__VERSION__ sim_method = sim_method spins = length(obj) time_points = length(t) adc_points=Ndims[1]
    @time timed_tuple = @timed run_sim_time_iter!(obj, seqd, sig, Xt, sim_method, itp; 
                                                  Nblocks, Nthreads, parts, w)
    # Result to CPU, if already in the CPU it does nothing
    sig = sum(sig; dims=length(Ndims)+1) |> cpu
    sig .*= get_adc_phase_compensation(seq)
    Xt = Xt |> cpu
    if enable_gpu GC.gc(true); CUDA.reclaim() end
    # Output
    if return_type == "state"
        out = Xt
    elseif return_type == "mat"
        out = sig
    elseif return_type == "raw"
        # To visually check the simulation blocks
        simParams_raw = copy(simParams)
        simParams_raw["sim_method"] = string(sim_method)
        simParams_raw["gpu"] = enable_gpu
        simParams_raw["Nthreads"] = Nthreads
        simParams_raw["t_sim_parts"] = t_sim_parts
        simParams_raw["Nblocks"] = Nblocks
        simParams_raw["sim_time_sec"] = timed_tuple.time
        out = signal_to_raw_data(sig, seq; phantom_name=obj.name, sys=sys, simParams=simParams_raw)
    end
    return out
end

"""
    M = simulate_slice_profile(seq; z, simParams)

Returns magnetization of spins distributed along `z` after running the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `z`: (`=range(-2e-2,2e-2,200)`) range for the z axis
- `simParams`: (`::Dict{String, Any}`, `=Dict{String,Any}("Δt_rf"=>1e-6)`) dictionary with
    simulation parameters

# Returns
- `M`: (`::Vector{Mag}`) final state of the Mag vector
"""
function simulate_slice_profile(seq; z=range(-2.e-2, 2.e-2, 200), simParams=Dict{String,Any}("Δt_rf" => 1e-6))
    simParams["return_type"] = "state"
    sys = Scanner()
    phantom = Phantom{Float64}(x=zeros(size(z)), z=Array(z))
    M = simulate(phantom, seq, sys; simParams)
    M
end
