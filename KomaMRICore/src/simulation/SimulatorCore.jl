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
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nthreads=Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads, type="ordered")
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        run_spin_precession!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
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
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nthreads=Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads; type="ordered")
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        run_spin_excitation!(@view(obj[p]), seq, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
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
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nblocks=1, Nthreads=Threads.nthreads(), parts=[1:length(seq)], excitation_bool=ones(Bool, size(parts)), w=nothing) where {T<:Real}
    # Simulation
    rfs = 0
    samples = 1
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        seq_block = @view seq[p]
        # Params
        # excitation_bool = is_RF_on(seq_block) #&& is_ADC_off(seq_block) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
        Nadc = sum(seq_block.ADC)
        acq_samples = samples:samples+Nadc-1
        dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times
        # Simulation wrappers
        if excitation_bool[block]
            run_spin_excitation_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
            rfs += 1
        else
            run_spin_precession_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
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

"""
Returns the result of a simple simulation using the disretization with blockvalues() function
"""
function komasim(seq::Sequence, obj::Phantom; Δtgr=1e-3, Δtrf=1e-5)

    # Create empty vectors to be filled during simulation for the sequence
    t, Δt, adct = Float64[], Float64[], Float64[]
    rfa, rfΔf, gxa, gya, gza = Float64[], Float64[], Float64[], Float64[], Float64[]
    rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask = Bool[], Bool[], Bool[], Bool[], Bool[]

    # Create empty vectors to be filled during simulation for the magnetizations and raw-signal
    magxy, magz, sig = Vector{ComplexF64}[], Vector{Float64}[], ComplexF64[]

    # Create the initial condition of the magnetization of the spins
    M_xy = zeros(ComplexF64, length(obj))
    M_z = obj.ρ

    # These are the initial times of the blocks
    to = cumsum([0.; durs(seq)])

    # Perform simulation iterating over every block-sequence
    for k = 1:length(seq)

        # Get the important vector values of the block-sequence
        tk, Δtk, rfak, rfΔfk, gxak, gyak, gzak, rf_onmaskk, gx_onmaskk, gy_onmaskk, gz_onmaskk, adc_onmaskk, adctk = blockvalues(seq, k; Δtgr, Δtrf)

        # Fill the vector values for the complete sequence
        # This first consideration is just for the initial condition of the complete sequence
        if k == 1
            push!(magxy, M_xy); push!(magz, M_z); push!(sig, sum(M_xy))
            append!(t, to[k] .+ tk[1]); append!(adct, to[k] .+ adctk)
            append!(rfa, rfak[1]); append!(rfΔf, rfΔfk[1]); append!(gxa, gxak[1]); append!(gya, gyak[1]); append!(gza, gzak[1])
            append!(rf_onmask, rf_onmaskk[1]); append!(adc_onmask, adc_onmaskk[1]); append!(gx_onmask, gx_onmaskk[1]); append!(gy_onmask, gy_onmaskk[1]); append!(gz_onmask, gz_onmaskk[1])
        end
        append!(t, to[k] .+ tk[2:end]); append!(Δt, Δtk); append!(adct, to[k] .+ adctk)
        append!(rfa, rfak[2:end]); append!(rfΔf, rfΔfk[2:end]); append!(gxa, gxak[2:end]); append!(gya, gyak[2:end]); append!(gza, gzak[2:end])
        append!(rf_onmask, rf_onmaskk[2:end]); append!(adc_onmask, adc_onmaskk[2:end]); append!(gx_onmask, gx_onmaskk[2:end]); append!(gy_onmask, gy_onmaskk[2:end]); append!(gz_onmask, gz_onmaskk[2:end])

        # Perform simulation iterating over each time step
        for i in eachindex(Δtk)

            # Excitation: Compute magnetization and signal when RF is on
            if rf_onmaskk[i]

                # B-field: compute the effective B field
                ΔBz = obj.Δw ./ (2π * γ) .- rfΔfk[i] ./ γ
                Bz = (gxak[i] .* obj.x .+ gyak[i] .* obj.y .+ gzak[i] .* obj.z) .+ ΔBz
                B = sqrt.((abs(rfak[i]))^2 .+ (abs.(Bz)).^2)
                B[B .== 0] .= eps(Float64)

                # Rotation: compute magnetization in rotation regime
                φ = (-2π * γ) * (B .* Δtk[i])
                nxy = rfak[i] ./ B
                nz = Bz ./ B
                α = cos.(φ/2) .- 1im*nz .* sin.(φ/2)
                β = -1im*nxy .* sin.(φ/2)
                Mxy = 2*conj.(α).*β.*M_z.+conj.(α).^2 .* M_xy.-β.^2 .*conj.(M_xy)
                Mz = (abs.(α).^2 .-abs.(β).^2).*M_z.-2*real.(α.*β.*conj.(M_xy))
                M_xy = Mxy
                M_z = Mz

                # Relaxation: compute magnetization in rotation regime
                M_xy = M_xy .* exp.(-Δtk[i] ./ obj.T2)
                M_z  = M_z  .* exp.(-Δtk[i] ./ obj.T1) .+ obj.ρ .* (1 .- exp.(-Δtk[i] ./ obj.T1))

            # Precession: compute magnetization and signal when RF is off
            else

                # B-field: compute effective B field
                ΔBz = obj.Δw ./ (2π * γ) .- rfΔfk[i] ./ γ
                Bz = (gxak[i] .* obj.x .+ gyak[i] .* obj.y .+ gzak[i] .* obj.z) .+ ΔBz
                ΔBzn = obj.Δw ./ (2π * γ) .- rfΔfk[i+1] ./ γ
                Bzn = (gxak[i+1] .* obj.x .+ gyak[i+1] .* obj.y .+ gzak[i+1] .* obj.z) .+ ΔBzn

                # Mxy: compute rotation and relaxation for Mxy in one step
                ϕ = (-2π * γ) * (0.5*(Bz+Bzn) .* Δtk[i])
                M_xy = M_xy .* exp.(1im .* (ϕ .- Δtk[i] ./ obj.T2))

                # Mz: compute just relaxation for Mz in one step (rotation phenomena doesn't happen)
                M_z  = M_z  .* exp.(-Δtk[i] ./ obj.T1) .+ obj.ρ .* (1 .- exp.(-Δtk[i] ./ obj.T1))

            end

            # Fill the magnetization and signal vectors
            push!(magxy, M_xy); push!(magz, M_z); push!(sig, sum(M_xy))

        end
    end

    return t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adct, reduce(hcat, magxy), reduce(hcat, magz), sig
end
