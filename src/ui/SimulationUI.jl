"""
    simulate!(w; sim_params=w.sim_params[])

Simulate the current window's phantom, sequence, scanner, and physiology. Update
`w.raw` and its view, save the reloadable MRD file, and return the raw data.
Blocks until complete; errors propagate and are recorded by `ui_status(w)`.
Rejects another load/simulation/reconstruction while busy. For asynchronous use,
keep a Julia task (`task = @async simulate!(w)`) and retrieve it with `fetch(task)`.
"""
function simulate!(w::KomaWindow; sim_params=w.sim_params[], initial=false)
    return ui_operation(w, :simulate) do
        run_simulation!(w, copy(sim_params); initial)
        return w.raw[]
    end
end

"""
    reconstruct!(w; rec_params=w.rec_params[])

Reconstruct the current window's raw data, update `w.img` and its view, and return
the reconstruction result. Blocks until complete; errors propagate and are recorded
by `ui_status(w)`. Rejects another load/simulation/reconstruction while busy.
Use `task = @async reconstruct!(w)` and `fetch(task)` for asynchronous callers.
Set `rec_params[:recon_policy]` to override label handling or coil combination,
for example `(COIL=:rss,)`.
"""
function reconstruct!(w::KomaWindow; rec_params=w.rec_params[], initial=false)
    return ui_operation(w, :reconstruct) do
        run_reconstruction!(w, copy(rec_params); initial)
        return w.img[]
    end
end

simulation_device_name(backend, selected) = selected
simulation_device_name(backend, ::Nothing) = KomaMRICore.device_name(backend)
simulation_device_name(backend, selected::Integer) = only(d.name for d in KomaMRICore._devices(backend) if d.id == selected)

function simulation_device_label(sim_params)
    use_gpu = Bool(get(sim_params, "gpu", true))
    backend = KomaMRICore.get_backend(use_gpu; verbose=false)
    backend_name = KomaMRICore.name(backend)
    if backend_name != "CPU"
        selected = get(sim_params, "gpu_device", nothing)
        device = simulation_device_name(backend, selected)
        return "GPU ($backend_name) · $device"
    end

    threads = get(sim_params, "Nthreads", Threads.nthreads())
    thread_label = "$threads thread$(threads == 1 ? "" : "s")"
    return "CPU · $(first(Sys.cpu_info()).model) · $thread_label"
end

function run_simulation!(w, sim_params; initial=false)
    message = initial ?
        "Precompiling and running simulation functions ..." : "Running simulation ..."
    params = KomaMRICore.default_sim_params(copy(sim_params))
    display_loading!(w, message)
    start_simulation_progress!(w)

    raw = try
        raw = fetch(Threads.@spawn begin
            backend = KomaMRICore.get_backend(params["gpu"]; verbose=false)
            if KomaMRICore.name(backend) != "CPU"
                isnothing(params["gpu_device"]) || KomaMRICore.set_device!(backend, params["gpu_device"])
                params["gpu_device"] = KomaMRICore.device_name(backend)
            end
            simulation_device = simulation_device_label(params)
            method = nameof(typeof(params["sim_method"]))
            status = "$simulation_device · $method · $(uppercase(params["precision"]))"
            display_loading!(w, message;
                details=DOM.div(DOM.div(status), DOM.br(),
                    DOM.div(DOM.div(), DOM.div(; class="fw-bold"); id="sim_timing")))
            simulate(
                w.obj[],
                w.seq[],
                w.sys[];
                sim_params,
                verbose=false,
                callbacks=(ui_progressbar_callback(w),),
                physio=w.physio[],
            )
        end)
        rawfile = joinpath(mktempdir(), "koma_sim.mrd")
        @info "Exporting to ISMRMRD file: $rawfile"
        save(ISMRMRDFile(rawfile), raw)
        w.files[:raw_data] = rawfile
        raw
    finally
        finish_simulation_progress!(w)
    end

    params = raw.params["userParameters"]
    sim_time = round(params["sim_time_sec"]; digits=3)
    backend_name = KomaMRICore.name(KomaMRICore.get_backend(params["gpu"]; verbose=false))
    simulation_backend = backend_name == "CPU" ? "CPU" : "GPU ($backend_name)"
    body = """
        <ul class="list-unstyled mb-0">
            <li><button type="button" class="btn koma-toast-action btn-circle btn-circle-sm m-1" title="View raw signal" aria-label="View raw signal" onclick="KomaUI.notify('view_raw_data')"><i class="bi bi-search"></i></button> Updating <b>Raw signal</b> plots ...</li>
            <li><button type="button" class="btn btn-primary btn-circle btn-circle-sm m-1" title="Reconstruct" aria-label="Reconstruct" onclick="KomaUI.notify('reconstruct')"><i class="bi bi-caret-right-fill"></i></button> Ready to <b>reconstruct</b>?</li>
        </ul>
    """
    update_filename!(w, "rawname", "koma_sim.mrd")
    toast!(w, 1, "$simulation_backend sim. successful<br>Time: $sim_time s", body)
    w.raw[] = raw
    return nothing
end

function run_reconstruction!(w, rec_params; initial=false)
    message = initial ?
        "Precompiling and running reconstruction functions ..." : "Running reconstruction ..."
    display_loading!(w, message)
    spinner = "<div class=\"spinner-border spinner-border-sm text-light\" role=\"status\"></div>"
    evaljs(w, js"""
        const button = document.getElementById('recon!');
        button.setAttribute('aria-busy', 'true');
        button.querySelector('.koma-action-label').innerHTML = $(spinner);
    """)

    try
        @info "Running reconstruction ..."
        reconstruction_result = fetch(Threads.@spawn @timed reconstruct_with_labels(w.raw[]; rec_params))
        body = """
            <ul class="list-unstyled mb-0"><li><button type="button" class="btn koma-toast-action btn-circle btn-circle-sm m-1" title="View reconstruction" aria-label="View reconstruction" onclick="KomaUI.notify('view_image')"><i class="bi bi-search"></i></button> Updating <b>Reconstruction</b> plots ...</li></ul>
        """
        rec_time = round(reconstruction_result.time; digits=3)
        warning = join(filter(!isempty, [reconstruction_geometry_warning(w.raw[]),
            reconstruction_geometry_warning(reconstruction_result.value)]), "\n")
        w.operation[] = UIOperation(; action=:reconstruct, state=:running, message=warning)
        isempty(warning) || (body = "<p style=\"white-space: pre-line;\">$(html_escape(warning))</p>" * body)
        toast!(w, 2, "Reconstruction successful<br>Time: $rec_time s", body;
            level=isempty(warning) ? :info : :warning)
        w.img[] = reconstruction_result.value
    finally
        evaljs(w, js"""
            const button = document.getElementById('recon!');
            button.querySelector('.koma-action-label').textContent = 'Reconstruct!';
            button.removeAttribute('aria-busy');
        """)
    end
    return nothing
end

function start_simulation_progress!(w::KomaWindow)
    isnothing(w.session[]) && return nothing
    evaljs(w, js"""(() => {
        const button = document.getElementById('simulate!');
        button.disabled = true;
        button.setAttribute('aria-busy', 'true');
        button.classList.add('koma-simulation-progress');
        button.style.setProperty('--simulation-progress', '0%');
        const bar = document.createElement('span');
        bar.id = 'simul_progress';
        bar.setAttribute('role', 'progressbar');
        bar.setAttribute('aria-label', 'Simulation progress');
        bar.setAttribute('aria-valuenow', '0');
        bar.setAttribute('aria-valuemin', '0');
        bar.setAttribute('aria-valuemax', '100');
        bar.textContent = '0%';
        button.querySelector('.koma-action-label').replaceChildren(bar);
        return true;
    })()""")
    return nothing
end

function finish_simulation_progress!(w::KomaWindow)
    evaljs(w, js"""
        const button = document.getElementById('simulate!');
        button.querySelector('.koma-action-label').textContent = 'Simulate!';
        button.removeAttribute('aria-busy');
        button.classList.remove('koma-simulation-progress');
        button.style.removeProperty('--simulation-progress');
        button.disabled = false;
    """)
    return nothing
end

function update_bonito_progress!(w::KomaWindow, block, Nblocks, elapsed_seconds, eta_seconds)
    progress = floor(Int, block / Nblocks * 100)
    w.operation[] = UIOperation(; action=:simulate, state=:running,
        progress=block / Nblocks, elapsed_seconds, eta_seconds)
    evaljs(w, js"""
        const bar = document.getElementById('simul_progress');
        if (bar) {
            bar.closest('button').style.setProperty('--simulation-progress', $progress + '%');
            bar.textContent = $progress + '%';
            bar.setAttribute('aria-valuenow', $progress);
        }
        const timing = document.getElementById('sim_timing');
        if (timing) {
            const formatDuration = seconds => {
                const total = Math.round(seconds);
                return Math.floor(total / 3600) + ' h ' +
                    Math.floor(total / 60) % 60 + ' min ' + total % 60 + ' s';
            };
            timing.children[0].textContent = 'Elapsed: ' + formatDuration($elapsed_seconds);
            timing.children[1].textContent = 'ETA: ' + formatDuration($eta_seconds);
        }
    """)
    return nothing
end

function ui_progressbar_callback(w::KomaWindow)
    started = time_ns()
    return Callback(1, (progress_info, _, _, _) -> begin
        elapsed_seconds = (time_ns() - started) / 1e9
        eta_seconds = elapsed_seconds * (progress_info.Nblocks - progress_info.block) / progress_info.block
        update_bonito_progress!(w, progress_info.block, progress_info.Nblocks, elapsed_seconds, eta_seconds)
    end)
end
