"""
    click!(w, action)

Activate a KomaUI menu or button action through the same handler as a UI click.
This runs the action directly, without finding or clicking a browser element.

Actions are `:home`, `:view_sequence`, `:view_kspace`, `:view_moment0`,
`:view_moment1`, `:view_moment2`, `:view_slew_rate`, `:view_phantom`, `:view_hardware_limits`,
`:view_coil_sensitivities`, `:view_raw_data`, `:view_image`, `:view_log_kspace`,
`:view_simulation_options`, `:view_reconstruction_options`, `:simulate`, and
`:reconstruct`. Reload actions are `:reload_sequence`, `:reload_phantom`,
`:reload_scanner`, and `:reload_raw_data`; they do nothing until a file is loaded.
MAT export actions are `:export_all`, `:export_sequence`, `:export_phantom`,
`:export_scanner`, `:export_raw_data`, and `:export_image`.
"""
function click!(w, action)
    event = string(action)
    w.handlers[event](event)
    return nothing
end

"""
    load_file!(w, target, filename)

Load a file as if selected in KomaUI's file picker. `target` is `:sequence`,
`:phantom`, `:scanner`, or `:raw_data` (explicit because `.h5` has multiple uses).
Updates the input, filename, toast, and plot, and remembers the path for reload.
Does not run simulation or reconstruction.
"""
function load_file!(w, target, filename)
    input = getproperty((sequence=w.seq, phantom=w.obj, scanner=w.sys, raw_data=w.raw), target)
    path = abspath(filename)
    value = callback_filepicker(path, w, input[])
    w.files[target] = path
    input[] = value
    return nothing
end

function handle(callback, w, event)
    w.handlers[event] = callback
    return nothing
end

function setup_actions!(w; darkmode)
    handle(w, "home") do _
        set_content!(w, w.home[], "index")
    end
    handle(w, "view_sequence") do _
        show_sequence!(w, w.seq[], :sequence; darkmode, physio=w.physio[])
    end
    handle(w, "view_kspace") do _
        show_sequence!(w, w.seq[], :kspace; darkmode)
    end
    handle(w, "view_moment0") do _
        show_sequence!(w, w.seq[], :moment0; darkmode)
    end
    handle(w, "view_moment1") do _
        show_sequence!(w, w.seq[], :moment1; darkmode)
    end
    handle(w, "view_moment2") do _
        show_sequence!(w, w.seq[], :moment2; darkmode)
    end
    handle(w, "view_slew_rate") do _
        show_sequence!(w, w.seq[], :slew_rate; darkmode, physio=w.physio[])
    end
    handle(w, "view_phantom") do _
        show_phantom!(w, w.obj[]; darkmode)
    end
    handle(w, "view_hardware_limits") do _
        show_scanner_parameters!(w, w.sys[])
    end
    handle(w, "view_coil_sensitivities") do _
        show_scanner!(w, w.sys[]; darkmode)
    end
    handle(w, "view_simulation_options") do _
        show_parameters!(w, w.sim_params[], "Simulation parameters", "simparams")
    end
    handle(w, "view_raw_data") do _
        show_signal!(w, w.raw[]; darkmode)
    end
    handle(w, "view_reconstruction_options") do _
        show_parameters!(w, w.rec_params[], "Reconstruction parameters", "recparams")
    end
    handle(w, "view_image") do _
        show_image!(w, w.img[], :absi; darkmode)
    end
    handle(w, "view_log_kspace") do _
        show_image!(w, w.img[], :absk; darkmode)
    end

    for target in (:sequence, :phantom, :scanner, :raw_data)
        handle(w, "reload_$target") do _
            haskey(w.files, target) && load_file!(w, target, w.files[target])
            return nothing
        end
    end
    for type in ("all", "sequence", "phantom", "scanner", "raw", "image")
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
        run_simulation!(w, w.sim_params[]; initial)
    end
    handle(w, "reconstruct") do _
        initial = is_first_rec
        is_first_rec = false
        run_reconstruction!(w, w.rec_params[]; initial)
    end
    return nothing
end
