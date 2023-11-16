abstract type SimulationMethod end #get all available types by using subtypes(KomaMRI.SimulationMethod)
abstract type SpinStateRepresentation{T<:Real} end #get all available types by using subtypes(KomaMRI.SpinStateRepresentation)

#Defined methods:
include("Bloch/BlochSimulationMethod.jl") #Defines Bloch simulation method
include("Bloch/BlochDictSimulationMethod.jl") #Defines BlochDict simulation method

function default_sim_params(simParams=Dict{String,Any}())
    get!(simParams, "gpu", true); if simParams["gpu"] check_use_cuda(); simParams["gpu"] &= use_cuda[] end
    get!(simParams, "gpu_device", 0)
    get!(simParams, "Nthreads", simParams["gpu"] ? 1 : Threads.nthreads())
    get!(simParams, "Nblocks", 20)
    get!(simParams, "Δt", 1e-3)
    get!(simParams, "Δt_rf", 5e-5)
    get!(simParams, "sim_method", Bloch())
    get!(simParams, "precision", "f32")
    get!(simParams, "return_type", "raw")
    simParams
end

# """
#     sig, Xt = run_spin_precession_parallel(obj, seq, M; Nthreads)

# Implementation in multiple threads for the simulation in free precession,
# separating the spins of the phantom `obj` in `Nthreads`.

# # Arguments
# - `obj`: (`::Phantom`) Phantom struct
# - `seq`: (`::Sequence`) Sequence struct

# # Keywords
# - `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
#     dividing the simulation into different phantom spin parts

# # Returns
# - `sig`: (`Vector{ComplexF64}`) raw signal over time
# - `Xt`: (`::Vector{Mag}`) final state of the Mag vector (or the initial state for the
#     next simulation step (the next step can be another precession step or an excitation
#     step))
# """
# function run_spin_precession_parallel!(obj::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
#     Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
#     Nthreads=Threads.nthreads()) where {T<:Real}

#     parts = kfoldperm(length(obj), Nthreads, type="ordered")
#     dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

#     ThreadsX.foreach(enumerate(parts)) do (i, p)
#         run_spin_precession!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
#     end

#     return nothing
# end

# """
#     M0 = run_spin_excitation_parallel(obj, seq, Xt; Nthreads)

# It gives rise to a rotation of M0 with an angle given by the efective magnetic field
# (including B1, gradients and off resonance) and with respect to a rotation axis. It uses
# different number threads to excecute the process.

# # Arguments
# - `obj`: (`::Phantom`) Phantom struct
# - `seq`: (`::Sequence`) Sequence struct

# # Keywords
# - `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
#     dividing the simulation into different phantom spin parts

# # Returns
# - `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (or the initial
#     state for the next precession simulation step)
# """
# function run_spin_excitation_parallel!(obj::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
#     Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
#     Nthreads=Threads.nthreads()) where {T<:Real}

#     parts = kfoldperm(length(obj), Nthreads; type="ordered")
#     dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

#     ThreadsX.foreach(enumerate(parts)) do (i, p)
#         run_spin_excitation!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
#     end

#     return nothing
# end

# """
#     S_interp, M0 = run_sim_time_iter(obj, seq, t, Δt; Nblocks, Nthreads, gpu, w)

# Performs the simulation over the total time vector `t` by dividing the time into `Nblocks`
# parts to reduce RAM usage and spliting the spins of the phantom `obj` into `Nthreads` to
# take advantage of CPU parallel processing.

# # Arguments
# - `obj`: (`::Phantom`) Phantom struct
# - `seq`: (`::Sequence`) Sequence struct
# - `t`: (`::Vector{Float64}`, `[s]`) non-uniform time vector
# - `Δt`: (`::Vector{Float64}`, `[s]`) delta time of `t`

# # Keywords
# - `Nblocks`: (`::Int`, `=16`) number of groups for spliting the simulation over time
# - `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
#     dividing the simulation into different phantom spin parts
# - `gpu`: (`::Function`) function that represents the gpu of the host
# - `w`: (`::Any`, `=nothing`) flag to regard a progress bar in the blink window UI. If
#     this variable is differnet from nothing, then the progress bar is considered

# # Returns
# - `S_interp`: (`::Vector{ComplexF64}`) interpolated raw signal
# - `M0`: (`::Vector{Mag}`) final state of the Mag vector
# """

# function run_sim_time_iter!(obj::Phantom, seq::DiscreteSequence, sig::AbstractArray{Complex{T}},
#     Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
#     Nblocks=1, Nthreads=Threads.nthreads(), parts=[1:length(seq)], excitation_bool=ones(Bool, size(parts)), w=nothing) where {T<:Real}
#     # Simulation
#     rfs = 0
#     samples = 1
#     progress_bar = Progress(Nblocks)
#     for (block, p) = enumerate(parts)
#         seq_block = @view seq[p]
#         # Params
#         # excitation_bool = is_RF_on(seq_block) #&& is_ADC_off(seq_block) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
#         Nadc = sum(seq_block.ADC)
#         acq_samples = samples:samples+Nadc-1
#         dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times
#         #println(dims)
#         #println(acq_samples)
#         # Simulation wrappers
#         if excitation_bool[block]
#             run_spin_excitation_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
#             rfs += 1
#         else
#             run_spin_precession_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
#         end
#         samples += Nadc
#         #Update progress
#         next!(progress_bar, showvalues=[(:simulated_blocks, block), (:rf_blocks, rfs), (:acq_samples, samples-1)])
#         update_blink_window_progress!(w, block, Nblocks)
#     end
#     return nothing
# end

# LUKE WORK

function run_sim_time_iter!(obj::Phantom, seq::DiscreteSequence, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nblocks=1, Nthreads=Threads.nthreads(), parts=[1:length(seq)], excitation_bool=ones(Bool, size(parts)), w=nothing) where {T<:Real}
    # Simulation
    rfs = 0
    samples = 1
    progress_bar = Progress(Nblocks)
    parts = kfoldperm(length(obj), Nthreads; type="ordered")
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times
    ThreadsX.foreach(enumerate(parts)) do (i, p)
        run_section_simulation!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
    end
    println("Simulation ended")
    return nothing
end

function run_section_simulation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    X = CuArray([1,0,0])
    Y = CuArray([0,1,0])
    Z = CuArray([0,0,1])
    points = xt' .* X .+ yt' .* Y .+ zt' .* Z
    M0 =  CUDA.real.(M.xy)' .* X .+ CUDA.imag.(M.xy)' .* Y.+ M.z' .* Z
    # Need to make B1 be a real complex vector (i.e mag and phase, not just mag)
    show(IOContext(stdout, :limit => true), "text/plain", M0)
    params = simParametersPhantom(Array(M.z)[1],p.T1,p.T2,1,γ,1.0,seq.Δt,points,M0,3.0,[0])
    bloch_symmetric_splitting!(real.(seq.B1), imag.(seq.B1),seq.Gx, seq.Gy, seq.Gz,params,sig,M,seq.ADC)
end   

# function run_section_simulation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
#     M::Mag{T}, sim_method::Bloch) where {T<:Real}
#     xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
#     yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
#     zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
#     points = xt' .* [1,0,0] .+ yt' .* [0,1,0] .+ zt' .* [0,0,1]
#     M0 =  real.(M.xy)' .* [1,0,0] .+ imag.(M.xy)' .* [0,1,0] .+ M.z' .* [0,0,1]
#     params = simParametersPhantom(1.0,p.T1,p.T2,1,γ,1.0,seq.Δt,points,M0,3.0,[0])
#     bloch_symmetric_splitting!(real.(seq.B1), imag.(seq.B1),seq.Gx, seq.Gy, seq.Gz,params,sig,seq.ADC)
# end   


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
    #Simulation parameter unpacking, and setting defaults if key is not defined
    simParams = default_sim_params(simParams)
    # Simulation init
    seqd = discretize(seq; simParams) # Sampling of Sequence waveforms
    parts, excitation_bool = get_sim_ranges(seqd; Nblocks=simParams["Nblocks"]) # Generating simulation blocks
    t_sim_parts = [seqd.t[p[1]] for p ∈ parts]; append!(t_sim_parts, seqd.t[end])
    # Spins' state init (Magnetization, EPG, etc.), could include modifications to obj (e.g. T2*)
    Xt, obj = initialize_spins_state(obj, simParams["sim_method"])
    # Signal init
    Ndims = sim_output_dim(obj, seq, sys, simParams["sim_method"])
    sig = zeros(ComplexF64, Ndims..., simParams["Nthreads"])
    # Objects to GPU
    if simParams["gpu"] #Default
        device!(simParams["gpu_device"])
        gpu_name = name.(devices())[simParams["gpu_device"]+1]
        obj  = obj  |> gpu #Phantom
        seqd = seqd |> gpu #DiscreteSequence
        Xt   = Xt   |> gpu #SpinStateRepresentation
        sig  = sig  |> gpu #Signal
    end
    if simParams["precision"] == "f32" #Default
        obj  = obj  |> f32 #Phantom
        seqd = seqd |> f32 #DiscreteSequence
        Xt   = Xt   |> f32 #SpinStateRepresentation
        sig  = sig  |> f32 #Signal
    elseif simParams["precision"] == "f64"
        obj  = obj  |> f64 #Phantom
        seqd = seqd |> f64 #DiscreteSequence
        Xt   = Xt   |> f64 #SpinStateRepresentation
        sig  = sig  |> f64 #Signal
    end
    # Simulation
    @info "Running simulation in the $(simParams["gpu"] ? "GPU ($gpu_name)" : "CPU with $(simParams["Nthreads"]) thread(s)")" koma_version=__VERSION__ sim_method = simParams["sim_method"] spins = length(obj) time_points = length(seqd.t) adc_points=Ndims[1]
    @time timed_tuple = @timed run_sim_time_iter!(obj, seqd, sig, Xt, simParams["sim_method"]; Nblocks=length(parts), Nthreads=simParams["Nthreads"], parts, excitation_bool, w)
    # Result to CPU, if already in the CPU it does nothing
    sig = sum(sig; dims=length(Ndims)+1) |> cpu
    sig .*= get_adc_phase_compensation(seq)
    Xt = Xt |> cpu
    if simParams["gpu"] GC.gc(true); CUDA.reclaim() end
    # Output
    if simParams["return_type"] == "state"
        out = Xt
    elseif simParams["return_type"] == "mat"
        out = sig
    elseif simParams["return_type"] == "raw"
        # To visually check the simulation blocks
        simParams_raw = copy(simParams)
        simParams_raw["sim_method"] = string(simParams["sim_method"])
        simParams_raw["gpu"] = simParams["gpu"]
        simParams_raw["Nthreads"] = simParams["Nthreads"]
        simParams_raw["t_sim_parts"] = t_sim_parts
        simParams_raw["type_sim_parts"] = excitation_bool
        simParams_raw["Nblocks"] = length(parts)
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
