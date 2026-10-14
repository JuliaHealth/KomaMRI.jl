"""
    ui_status(w)

Return a snapshot of the latest load, simulation, reconstruction, or device selection:
`action`, `state` (`:idle`, `:running`, `:succeeded`, or `:failed`), `progress`
(fraction in `[0, 1]`, or `nothing` when unknown), and `message` (warning or error text).
`elapsed_seconds` and `eta_seconds` are simulation-only (`nothing` for other
operations or before the first update). ETA is a block-based estimate and can
fluctuate during warmup.
Also includes the current `view`, a copy of loaded `files`, browser `url` (or `nothing`),
and whether a browser/desktop session is `connected`.
Completion means Julia data and view updates have finished; it does not wait for
the browser to finish painting the plot.
"""
function ui_status(w)
    operation = w.operation[]
    return (; operation.action, operation.state, operation.progress,
        operation.elapsed_seconds, operation.eta_seconds, operation.message,
        view=w.state[], files=copy(w.files),
        url=isnothing(w.server[]) ? nothing : Bonito.HTTPServer.online_url(w.server[], window_route(w)),
        connected=!isnothing(w.session[]) && isopen(w.session[]))
end

"""
    ui_actions(w)

Return the supported `click!` action names as a sorted vector of symbols.
"""
ui_actions(w) = sort!(Symbol.(collect(keys(w.handlers))))

function ui_operation(callback, w, action)
    trylock(w.operation_lock) || error("KomaUI is busy with $(w.operation[].action).")
    if w.operation[].state === :running
        unlock(w.operation_lock)
        error("KomaUI is busy with $(w.operation[].action).")
    end
    previous_content, previous_state = w.content[], w.state[]
    try
        w.operation[] = UIOperation(; action, state=:running)
        result = callback()
        operation = w.operation[]
        w.operation[] = UIOperation(; action, state=:succeeded, progress=1.0,
            elapsed_seconds=operation.elapsed_seconds, eta_seconds=operation.eta_seconds,
            message=operation.message)
        return result
    catch error
        w.operation[] = UIOperation(; action, state=:failed,
            elapsed_seconds=w.operation[].elapsed_seconds, message=sprint(showerror, error))
        restore_content!(w, previous_content, previous_state)
        rethrow()
    finally
        unlock(w.operation_lock)
    end
end

function ui_action(callback, w, action)
    try
        callback()
    catch error
        @error "KomaUI action failed" action exception=(error, catch_backtrace())
        failure_toast!(w, 0, replace(string(action), '_' => ' '), error)
    end
    return nothing
end

"""
    click!(w, action; kwargs...)

Activate a KomaUI menu or button action through the same handler as a UI click.
This runs the action directly, without finding or clicking a browser element.
It is synchronous and propagates errors. Use `simulate!`, `reconstruct!`, or
`load_file!` when the returned data is needed; `click!` returns `nothing`.
Use `ui_actions(w)` to discover actions and `ui_status(w)` to inspect operation status.

Actions are `:home`, `:toggle_sidebar`, `:toggle_theme`, `:view_sequence`, `:view_sequence_definitions`, `:view_kspace`, `:view_moment0`,
`:view_moment1`, `:view_moment2`, `:view_slew_rate`, `:view_phantom`, `:view_hardware_limits`,
`:view_phantom_slice`, `:view_coil_sensitivities`, `:view_raw_data`, `:view_image`, `:view_log_kspace`,
`:view_configuration`, `:view_simulation_options`, `:view_reconstruction_options`, `:simulate`, and
`:reconstruct`. Reload actions are `:reload_sequence`, `:reload_phantom`,
`:reload_scanner`, `:reload_raw_data`, and `:reload_configuration`; they do nothing until a file is loaded.
MAT export actions are `:export_all`, `:export_sequence`, `:export_phantom`,
`:export_scanner`, `:export_raw_data`, and `:export_image`.
Export raw data as ISMRMRD with `:export_mrd`.
Use `export_results(w, output; path=...)` to save without an interactive picker.

Select the next simulation's device with `click!(w, :select_device; gpu=false)`
for CPU or `click!(w, :select_device; gpu=true, id=device_id)` for a GPU listed by
`get_devices()`. Omitting `id` selects the backend's current GPU. GPU selection
requires a loaded, functional backend; unavailable devices throw without changing
the settings. Selection is recorded by `ui_status(w)` and rejected while busy.
"""
function click!(w, action; kwargs...)
    event = string(action)
    if w.operation[].state === :running && !(event in ("toggle_theme", "toggle_sidebar"))
        error("KomaUI is busy with $(w.operation[].action).")
    end
    w.handlers[event](event; kwargs...)
    return nothing
end

"""
    load_file!(w, target, filename)

Load a file as if selected in KomaUI's file picker. `target` is `:sequence`,
`:phantom`, `:scanner`, `:raw_data`, or `:configuration` (explicit because `.h5` has multiple uses).
Configuration TOML files merge `[sim_params]` and `[recon_params]` into the window's
current parameters; omitted settings are unchanged. Both sections are parsed and
converted before applying either. Returns the merged parameter dictionaries.
Updates the input, filename, toast, and plot, and remembers the path for reload.
Returns the loaded object after updating the view and checking sequence geometry
and ADC labels. Warnings are displayed in the UI; loading/check failures throw.
Calls block until complete and reject concurrent loads, simulation, or reconstruction.
Does not run simulation or reconstruction.
"""
function load_file!(w, target, filename)
    return ui_operation(w, Symbol("load_", target)) do
        path = abspath(filename)
        if target === :configuration
            params = read_configuration(path; sim_params=w.sim_params[], recon_params=w.rec_params[])
            params.sim_params == w.sim_params[] || (w.sim_params[] = params.sim_params)
            params.recon_params == w.rec_params[] || (w.rec_params[] = params.recon_params)
            show_configuration!(w)
            w.files[target] = path
            update_filename!(w, "configname", basename(path))
            toast!(w, 0, "Loaded configuration", html_escape(basename(path)))
            return params
        end
        input = getproperty((sequence=w.seq, phantom=w.obj, scanner=w.sys, raw_data=w.raw), target)
        value = callback_filepicker(path, w, input[])
        w.files[target] = path
        input[] = value
        if target === :sequence
            warning = fetch(Threads.@spawn sequence_reconstruction_warning(value))
            w.operation[] = UIOperation(; action=Symbol("load_", target), state=:running, message=warning)
            isempty(warning) || loaded_data_toast(w, basename(path), "view_sequence", "simulate", "Sequence"; warning)
        end
        return value
    end
end

function handle(callback, w, event)
    w.handlers[event] = callback
    return nothing
end

function setup_actions!(w)
    handle(function (event; gpu, id=nothing)
        ui_operation(w, :select_device) do
            devices = filter(d -> (d.backend != "CPU") == gpu, get_devices(gpu))
            isempty(devices) && throw(ArgumentError("No functional GPU backend is loaded."))
            index = findfirst(d -> isnothing(id) ? d.current : d.id == id, devices)
            isnothing(index) && throw(ArgumentError("Device $id is unavailable. Choose an ID from get_devices()."))
            device = devices[index]
            w.sim_params[] = merge(w.sim_params[], Dict("gpu" => gpu, "gpu_device" => device.id))
        end
        return nothing
    end, w, "select_device")
    handle(w, "home") do _
        set_content!(w, w.home[], "index")
    end
    handle(w, "toggle_sidebar") do _
        evaljs(w, js"toggleSidebar()")
    end
    handle(w, "toggle_theme") do _
        w.darkmode[] = !w.darkmode[]
    end
    handle(w, "view_sequence_definitions") do _
        show_parameters!(w, w.seq[].DEF, "Sequence definitions", "seqdefs")
    end
    for view in (:sequence, :kspace, :moment0, :moment1, :moment2, :slew_rate)
        handle(w, "view_$view") do _
            show_sequence!(w, w.seq[], view; physio=w.physio[])
        end
    end
    handle(w, "view_phantom") do _
        show_phantom!(w, w.obj[])
    end
    handle(w, "view_phantom_slice") do _
        show_phantom!(w, w.obj[]; slice_view=true)
    end
    handle(w, "view_hardware_limits") do _
        show_scanner_parameters!(w, w.sys[])
    end
    handle(w, "view_coil_sensitivities") do _
        show_scanner!(w, w.sys[])
    end
    handle(w, "view_simulation_options") do _
        show_parameters!(w, w.sim_params[], "Simulation parameters", "simparams")
    end
    handle(w, "view_configuration") do _
        show_configuration!(w)
    end
    handle(w, "view_raw_data") do _
        show_signal!(w, w.raw[])
    end
    handle(w, "view_reconstruction_options") do _
        show_parameters!(w, w.rec_params[], "Reconstruction parameters", "recparams")
    end
    handle(w, "view_image") do _
        show_image!(w, w.img[], :absi)
    end
    handle(w, "view_log_kspace") do _
        show_image!(w, w.img[], :absk)
    end

    for target in (:sequence, :phantom, :scanner, :raw_data, :configuration)
        handle(w, "reload_$target") do _
            haskey(w.files, target) && load_file!(w, target, w.files[target])
            return nothing
        end
    end
    for type in ("all", "sequence", "phantom", "scanner", "raw", "image", "mrd")
        event = type == "raw" ? "export_raw_data" : "export_$type"
        handle(w, event) do _
            save_ui!(w, w.seq[], w.obj[], w.sys[], w.raw[], w.img[], w.rec_params[]; type)
        end
    end

    is_first_sim = true
    is_first_rec = true
    handle(w, "simulate") do _
        initial = is_first_sim
        is_first_sim = false
        simulate!(w; initial)
    end
    handle(w, "reconstruct") do _
        initial = is_first_rec
        is_first_rec = false
        reconstruct!(w; initial)
    end
    return nothing
end
