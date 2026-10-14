using PrecompileTools: @setup_workload, @compile_workload
using Preferences: load_preference
import TOML

const CLI_BACKENDS = ("AMDGPU", "CUDA", "Metal", "oneAPI")

const CLI_ALIASES = Dict(
    "-i" => "--inputs",
    "-o" => "--outputs",
    "-b" => "--backend",
    "-s" => "--sim-param",
    "-r" => "--recon-param",
    "-h" => "--help",
)

const CLI_HELP = """
KomaMRI command line app.

Usage:
  koma
  koma --browser --port 9384
  koma -i epi.seq brain.phantom
  koma -i epi.seq brain.phantom -o raw.mrd [image.mat]

Options:
  -i, --inputs FILE...             .seq, .phantom/.h5, .sys in any order
  -o, --outputs RAW|_ [IMAGE.mat]  RAW is .mrd or .mat; use _ for recon-only
  -b, --backend NAME               CPU (default), CUDA, Metal, AMDGPU, or oneAPI
  -s, --sim-param KEY=VALUE        repeatable simulation parameter
  -r, --recon-param KEY=VALUE      repeatable reconstruction parameter
      --config FILE               TOML simulation/reconstruction parameters
      --browser                    serve the UI without opening Electron
      --host ADDRESS               browser bind address (default: 127.0.0.1)
      --port PORT                  browser server port (default: 9384)
  -h, --help                       show this help

Julia flags:
  koma --threads=8 -- [koma args...]
"""

Base.@kwdef mutable struct CLIOptions
    sequence::Union{Nothing,String} = nothing
    phantom::Union{Nothing,String} = nothing
    scanner::Union{Nothing,String} = nothing
    sim_output::Union{Nothing,String} = nothing
    recon_output::Union{Nothing,String} = nothing
    backend::String = "CPU"
    browser::Bool = false
    host::String = "127.0.0.1"
    port::Int = 9384
    sim_params::Dict{String,Any} = Dict{String,Any}("gpu" => false)
    recon_params::Dict{Symbol,Any} = Dict{Symbol,Any}(:reco => "direct")
end

function parse_cli_value(x::AbstractString)
    x in ("true", "false") && return x == "true"
    y = tryparse(Int, x)
    isnothing(y) || return y
    z = tryparse(Float64, x)
    isnothing(z) || return z
    if startswith(x, "[") && endswith(x, "]")
        return parse_cli_value.(strip.(split(x[2:end-1], ",")))
    end
    if startswith(x, "(") && endswith(x, ")")
        return Tuple(parse_cli_value.(strip.(split(x[2:end-1], ","))))
    end
    return String(x)
end

function parse_cli_sim_method(x::AbstractString)
    sym = Symbol(x)
    isdefined(KomaMRICore, sym) || error("Unsupported simulation method: $x")
    sim_method = getproperty(KomaMRICore, sym)
    if !(sim_method isa Type && sim_method <: KomaMRICore.SimulationMethod)
        error("Unsupported simulation method: $x")
    end
    return sim_method()
end

function parse_cli_value(key::AbstractString, x::AbstractString)
    key == "sim_method" && return parse_cli_sim_method(x)
    return parse_cli_value(x)
end

param_key(::Dict{String,Any}, key) = String(key)
param_key(::Dict{Symbol,Any}, key) = Symbol(key)

function parse_cli_param!(params, param::AbstractString)
    key, value = split(param, "="; limit=2)
    params[param_key(params, key)] = parse_cli_value(key, value)
    return params
end

function cli_param_value(key, value)
    key == "sim_method" && value isa AbstractString && return parse_cli_sim_method(value)
    key == "recon_policy" && return (; (Symbol(k) => Symbol(v) for (k, v) in pairs(value))...)
    return value
end

function read_configuration(filename; sim_params=Dict{String,Any}(), recon_params=Dict{Symbol,Any}())
    config = TOML.parsefile(filename)
    all(in(("sim_params", "recon_params")), keys(config)) ||
        throw(ArgumentError("Configuration accepts [sim_params] and [recon_params] sections."))
    sim = Dict{String,Any}(key => cli_param_value(key, value)
        for (key, value) in get(config, "sim_params", Dict()))
    rec = Dict{Symbol,Any}(Symbol(key) => cli_param_value(key, value)
        for (key, value) in get(config, "recon_params", Dict()))
    if haskey(rec, :recon_policy)
        rec[:recon_policy] = reconstruction_policy(merge(
            get(recon_params, :recon_policy, NamedTuple()), rec[:recon_policy]))
    end
    return (; sim_params=merge(sim_params, sim), recon_params=merge(recon_params, rec))
end

function merge_cli_preferences!(opts, prefs)
    haskey(prefs, "backend") && (opts.backend = String(prefs["backend"]))

    inputs = get(prefs, "inputs", Dict{String,Any}())
    haskey(inputs, "sequence") && (opts.sequence = String(inputs["sequence"]))
    haskey(inputs, "phantom") && (opts.phantom = String(inputs["phantom"]))
    haskey(inputs, "scanner") && (opts.scanner = String(inputs["scanner"]))

    outputs = get(prefs, "outputs", Dict{String,Any}())
    haskey(outputs, "rawdata") && (opts.sim_output = String(outputs["rawdata"]))
    haskey(outputs, "image") && (opts.recon_output = String(outputs["image"]))

    for (key, value) in get(prefs, "sim_params", Dict{String,Any}())
        key = String(key)
        opts.sim_params[key] = cli_param_value(key, value)
    end
    for (key, value) in get(prefs, "recon_params", Dict{String,Any}())
        key = String(key)
        opts.recon_params[Symbol(key)] = cli_param_value(key, value)
    end
    return opts
end

function load_cli_preferences!(opts)
    return merge_cli_preferences!(opts, load_preference(@__MODULE__, "koma", Dict{String,Any}()))
end

is_cli_flag(arg) = startswith(arg, "-")

function split_cli_option(arg)
    parts = split(arg, "="; limit=2)
    name = get(CLI_ALIASES, parts[1], parts[1])
    value = length(parts) == 2 ? String(parts[2]) : nothing
    return name, value
end

function take_cli_value(args, i, name, value)
    isnothing(value) || return value, i + 1
    i == length(args) && error("Missing value for $name")
    return args[i + 1], i + 2
end

function take_cli_values(args, i, name, value; max)
    isnothing(value) || error("$name does not accept =")
    values = String[]
    i += 1
    while i <= length(args) && !is_cli_flag(args[i])
        push!(values, args[i])
        i += 1
    end
    isempty(values) && error("Missing value for $name")
    length(values) <= max || error("$name accepts at most $max values")
    return values, i
end

function parse_cli_input!(opts, input)
    input == "_" && return opts
    ext = splitext(input)[2]
    if ext == ".seq"
        opts.sequence = input
    elseif ext in (".phantom", ".h5")
        opts.phantom = input
    elseif ext == ".sys"
        opts.scanner = input
    else
        error("Unsupported input extension: $ext")
    end
    return opts
end

function parse_cli_outputs!(opts, outputs)
    outputs == ["_"] && return opts
    if length(outputs) == 1
        opts.sim_output = only(outputs)
    elseif outputs[1] == "_"
        opts.recon_output = outputs[2]
    elseif splitext(outputs[2])[2] == ".mrd"
        opts.sim_output = outputs[2]
        outputs[1] == "_" || (opts.recon_output = outputs[1])
    else
        opts.sim_output = outputs[1]
        outputs[2] == "_" || (opts.recon_output = outputs[2])
    end
    return opts
end

function parse_cli_args(args, opts=CLIOptions())
    configuration_file = nothing
    sim_overrides, recon_overrides = Dict{String,Any}(), Dict{Symbol,Any}()
    i = 1
    while i <= length(args)
        name, value = split_cli_option(args[i])
        if name == "--inputs"
            inputs, i = take_cli_values(args, i, name, value; max=3)
            foreach(input -> parse_cli_input!(opts, input), inputs)
        elseif name == "--outputs"
            outputs, i = take_cli_values(args, i, name, value; max=2)
            parse_cli_outputs!(opts, outputs)
        elseif name == "--backend"
            opts.backend, i = take_cli_value(args, i, name, value)
        elseif name == "--sim-param"
            param, i = take_cli_value(args, i, name, value)
            parse_cli_param!(sim_overrides, param)
        elseif name == "--recon-param"
            param, i = take_cli_value(args, i, name, value)
            parse_cli_param!(recon_overrides, param)
        elseif name == "--config"
            configuration_file, i = take_cli_value(args, i, name, value)
        elseif name == "--browser"
            opts.browser = isnothing(value) ? true : parse(Bool, value)
            i += 1
        elseif name == "--host"
            opts.host, i = take_cli_value(args, i, name, value)
        elseif name == "--port"
            port, i = take_cli_value(args, i, name, value)
            opts.port = parse(Int, port)
        else
            error("Unknown argument: $(args[i])")
        end
    end
    if !isnothing(configuration_file)
        params = read_configuration(configuration_file; sim_params=opts.sim_params, recon_params=opts.recon_params)
        opts.sim_params, opts.recon_params = params.sim_params, params.recon_params
    end
    merge!(opts.sim_params, sim_overrides)
    merge!(opts.recon_params, recon_overrides)
    return opts
end

function load_cli_backend!(opts)
    if opts.backend == "CPU"
        opts.sim_params["gpu"] = false
    elseif opts.backend in CLI_BACKENDS
        # App shims pin LOAD_PATH; expose the default env for optional GPU backends.
        "@v#.#" in LOAD_PATH || push!(LOAD_PATH, "@v#.#")
        Base.eval(Main, :(using $(Symbol(opts.backend))))
        opts.sim_params["gpu"] = true
    else
        error("Unsupported backend: $(opts.backend)")
    end
    return nothing
end

function cli_inputs(opts)
    sys = isnothing(opts.scanner) ? setup_scanner() : read_scanner(opts.scanner)
    seq = isnothing(opts.sequence) ? setup_sequence(sys) : read_seq(opts.sequence)
    obj = isnothing(opts.phantom) ? setup_phantom() : load_cli_phantom(opts.phantom)
    return sys, seq, obj
end

function load_cli_phantom(filename)
    ext = splitext(filename)[2]
    ext == ".phantom" && return read_phantom(filename)
    ext == ".h5" && return read_phantom_jemris(filename)
    error("Unsupported phantom extension: $ext")
end

function cli_output_dir(filename)
    dir = dirname(filename)
    return isempty(dir) ? "." : dir
end

function mk_cli_output_dir(filename)
    dir = cli_output_dir(filename)
    mkpath(dir)
    return dir
end

function save_cli_raw(raw, filename)
    dir = mk_cli_output_dir(filename)
    ext = splitext(filename)[2]
    if ext == ".mrd"
        save(ISMRMRDFile(filename), raw)
    elseif ext == ".mat"
        export_2_mat_raw(raw, dir; matfilename=basename(filename))
    else
        error("Unsupported simulation output extension: $ext")
    end
    return nothing
end

function reconstruct_cli(raw, rec_params)
    return reconstruct_with_labels(raw; rec_params)
end

function save_cli_recon(image, rec_params, filename)
    splitext(filename)[2] == ".mat" || error("Unsupported reconstruction output extension: $(splitext(filename)[2])")
    dir = mk_cli_output_dir(filename)
    export_2_mat_image(image, rec_params, dir; matfilename=basename(filename))
    return nothing
end

function print_cli_versions()
    @info "KomaMRI loaded successfully 🚀" KomaMRI=string(pkgversion(@__MODULE__)) KomaMRIBase=string(pkgversion(KomaMRIBase)) KomaMRICore=string(pkgversion(KomaMRICore)) KomaMRIFiles=string(pkgversion(KomaMRIFiles)) KomaMRIPlots=string(pkgversion(KomaMRIPlots))
    return nothing
end

function keep_app_open(w)
    wait(w)
    return nothing
end

function run_cli(args)
    opts = parse_cli_args(args, load_cli_preferences!(CLIOptions()))
    load_cli_backend!(opts)
    # Loading a GPU backend may activate package extensions in a newer world.
    return Base.invokelatest(run_cli, opts)
end

function run_cli(opts::CLIOptions)
    sys, seq, obj = cli_inputs(opts)
    if isnothing(opts.sim_output) && isnothing(opts.recon_output)
        w = KomaUI(; sys, seq, obj, sim=opts.sim_params, rec=opts.recon_params,
            browser=opts.browser, host=opts.host, port=opts.port, verbose=false, return_window=true)
        keep_app_open(w)
        return nothing
    end

    print_cli_versions()
    raw = simulate(obj, seq, sys; sim_params=opts.sim_params)
    isnothing(opts.sim_output) || save_cli_raw(raw, opts.sim_output)
    if !isnothing(opts.recon_output)
        image = reconstruct_cli(raw, opts.recon_params)
        save_cli_recon(image, opts.recon_params, opts.recon_output)
    end
    return nothing
end

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            opts = CLIOptions()
            sys, seq, obj = cli_inputs(opts)
            w = KomaUI(;
                sys,
                seq,
                obj,
                sim=opts.sim_params,
                rec=opts.recon_params,
                show_window=false,
                return_window=true,
                verbose=false,
            )
            close(w)
        end
    end
end
