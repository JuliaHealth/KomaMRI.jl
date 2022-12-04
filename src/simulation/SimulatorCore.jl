const Nphyscores = Hwloc.num_physical_cores()

abstract type SimulationMethod end #get all available types by using subtypes(KomaMRI.SimulationMethod)
abstract type SpinStateRepresentation{T<:Real} end #get all available types by using subtypes(KomaMRI.SpinStateRepresentation)

#Defined methods:
include("Bloch/BlochSimulationMethod.jl") #Defines Bloch simulation method

"""
    sig, Xt = run_spin_precession_parallel(obj, seq, M; Nthreads)

Implementation in multiple threads for the simulation in free precession,
separating the spins of the phantom `obj` in `Nthreads`.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `sig`: (`Vector{ComplexF64}`) raw signal over time
- `Xt`: (`::Vector{Mag}`) final state of the Mag vector (or the initial state for the
    next simulation step (the next step can be another precession step or an excitation
    step))
"""
NVTX.@range function run_spin_precession_parallel(obj::Phantom{T}, seq::DiscreteSequence{T}, Xt::SpinStateRepresentation{T};
    Nthreads=Nphyscores) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads, type="ordered")

    sig = ThreadsX.mapreduce(+, parts) do p #Thread-safe summation
        sig_p, Xt[p] = run_spin_precession(obj[p], seq, Xt[p])
        sig_p
    end

    return sig, Xt
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
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (or the initial
    state for the next precession simulation step)
"""
NVTX.@range function run_spin_excitation_parallel(obj::Phantom{T}, seq::DiscreteSequence{T}, Xt::SpinStateRepresentation{T};
    Nthreads=Nphyscores) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads; type="ordered")

    Threads.@threads for p ∈ parts
        Xt[p] = run_spin_excitation(obj[p], seq, Xt[p])
    end

    return Xt
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
- `Nthreads`: (`::Int`, `=Hwloc.num_physical_cores()`) number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) function that represents the gpu of the host
- `w`: (`::Any`, `=nothing`) flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `S_interp`: (`::Vector{ComplexF64}`) interpolated raw signal
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
NVTX.@range function run_sim_time_iter(obj::Phantom, seq::DiscreteSequence, Xt::SpinStateRepresentation{T}, sig::AbstractArray{Complex{T}};
    Nblocks=1, Nthreads=Nphyscores, parts=[1:length(seq)], w=nothing) where {T<:Real}
    # Simulation
    rfs = 0
    samples = 1
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        # params
        excitation_bool = is_RF_on(seq[p]) && is_ADC_off(seq[p]) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
        Nadc = sum(seq[p].ADC)
        # Simulation wrappers
        if excitation_bool
            Xt = run_spin_excitation_parallel(obj, seq[p], Xt; Nthreads)
            rfs += 1
        else
            sig_aux, Xt = run_spin_precession_parallel(obj, seq[p], Xt; Nthreads)
            sig[samples:samples+Nadc-1, :] .= sig_aux'
            samples += Nadc
        end
        #Update progress
        next!(progress_bar, showvalues=[(:simulated_blocks, block), (:rf_blocks, rfs), (:adc_samples, samples)])
        update_blink_window_progress!(w, block, Nblocks)
    end
    return sig
end

"""Updates KomaUI's simulation progress bar."""
function update_blink_window_progress!(w, block, Nblocks)
    if w !== nothing #update Progress to Blink Window
        progress = string(floor(Int, block / Nblocks * 100))
        @js_ w (@var progress = $progress;
        document.getElementById("simul_progress").style.width = progress + "%";
        document.getElementById("simul_progress").innerHTML = progress + "%";
        document.getElementById("simul_progress").setAttribute("aria-valuenow", progress))
    end
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
julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq("examples/1.sequences/spiral.seq");

julia> ismrmrd = simulate(obj, seq, sys);

julia> plot_signal(ismrmrd)
```
"""
NVTX.@range function simulate(obj::Phantom, seq::Sequence, sys::Scanner; simParams=Dict{String,Any}(), w=nothing)
    #Simulation parameter parsing
    enable_gpu  = get(simParams, "gpu", use_cuda)
    Nthreads    = get(simParams, "Nthreads", enable_gpu ? 1 : Nphyscores)
    Δt          = get(simParams, "Δt", 1e-3)
    Δt_rf       = get(simParams, "Δt_rf", 1e-4)
    return_type = get(simParams, "return_type", "raw")
    Nblocks     = get(simParams, "Nblocks", 1)
    sim_method  = get(simParams, "sim_method", Bloch())
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
    ADCflag = [any(tt .== tadc) for tt ∈ t[2:end]]
    seqd = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    # Spins' state init (Magnetization, EPG, etc.)
    Xt = initialize_spins_state(obj, sim_method)
    # Signal init
    Nadc = sum(seq.ADC.N)
    Ncoils = 1
    sig = zeros(ComplexF64, Nadc, Ncoils)
    # Objects to GPU
    if enable_gpu
        d = 0
        device!(d)
        gpu_name = name.(devices())[d+1]
        obj  = obj  |> gpu #Phantom
        seqd = seqd |> gpu #DiscreteSequence
        Xt   = Xt   |> gpu #SpinStateRepresentation
        sig  = sig  |> gpu #Signal
    end
    # Simulation
    @info "Running simulation in the $(enable_gpu ? "GPU ($gpu_name)" : "CPU with $Nthreads thread(s)")" sim_method = sim_method spins = length(obj) time_points = length(t) adc_points=Nadc
    timed_tuple = @timed run_sim_time_iter(obj, seqd, Xt, sig; Nblocks, Nthreads, parts, w)
    @info "Simulation finished in $(timed_tuple.time) seconds!"
    #Result to CPU, if already in the CPU it does nothing
    sig = timed_tuple.value |> cpu
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
        simParams_raw["sim_time"] = timed_tuple.time
        out = signal_to_raw_data(sig, seq; phantom=obj, sys=sys, simParams=simParams_raw)
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
function simulate_slice_profile(seq; z=range(-2e-2, 2e-2, 200), simParams=Dict{String,Any}("Δt_rf" => 1e-6))
    simParams["return_type"] = "state"
    sys = Scanner()
    phantom = Phantom(; x=zeros(size(z)), z)
    M = simulate(phantom, seq, sys; simParams)
    M
end
