function simulation_device_label(sim_params)
    backend = KomaMRICore.get_backend(Bool(get(sim_params, "gpu", true)))
    backend_name = KomaMRICore.name(backend)
    backend_name == "CPU" || return "GPU ($backend_name)"

    threads = get(sim_params, "Nthreads", Threads.nthreads())
    return "CPU ($threads thread$(threads == 1 ? "" : "s"))"
end

function run_simulation!(w, sim_params; initial=false)
    previous_content = w.content[]
    previous_state = w.state[]
    message = initial ?
        "Precompiling and running simulation functions ..." : "Running simulation ..."
    params = KomaMRICore.default_sim_params(copy(sim_params))
    simulation_device = simulation_device_label(params)
    method = nameof(typeof(params["sim_method"]))
    status = "$simulation_device · $method · $(uppercase(params["precision"]))"
    display_loading!(w, message; details=status)
    start_simulation_progress!(w)

    raw = try
        raw = simulate(
            w.obj[],
            w.seq[],
            w.sys[];
            sim_params,
            callbacks=(ui_progressbar_callback(w),),
            physio=w.physio[],
        )
        rawfile = joinpath(mktempdir(), "koma_sim.mrd")
        @info "Exporting to ISMRMRD file: $rawfile"
        save(ISMRMRDFile(rawfile), raw)
        w.files[:raw_data] = rawfile
        raw
    catch error
        @error "Simulation failed" exception=(error, catch_backtrace())
        restore_content!(w, previous_content, previous_state)
        failure_toast!(w, 1, "Simulation", error)
        return nothing
    finally
        finish_simulation_progress!(w)
    end

    params = raw.params["userParameters"]
    sim_time = round(params["sim_time_sec"]; digits=3)
    simulation_device = simulation_device_label(params)
    body = """
        <ul class="list-unstyled mb-0">
            <li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View raw signal" aria-label="View raw signal" onclick="KomaUI.notify('view_raw_data')"><i class="bi bi-search"></i></button> Updating <b>Raw signal</b> plots ...</li>
            <li><button type="button" class="btn btn-primary btn-circle btn-circle-sm m-1" title="Reconstruct" aria-label="Reconstruct" onclick="KomaUI.notify('reconstruct')"><i class="bi bi-caret-right-fill"></i></button> Ready to <b>reconstruct</b>?</li>
        </ul>
    """
    update_filename!(w, "rawname", "koma_sim.mrd")
    toast!(w, 1, "$simulation_device sim. successful<br>Time: $sim_time s", body)
    w.raw[] = raw
    return nothing
end

function run_reconstruction!(w, rec_params; initial=false)
    previous_content = w.content[]
    previous_state = w.state[]
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
        reconstruction_result = @timed reconstruct_with_labels(w.raw[]; rec_params)
        body = """
            <ul class="list-unstyled mb-0"><li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View reconstruction" aria-label="View reconstruction" onclick="KomaUI.notify('view_image')"><i class="bi bi-search"></i></button> Updating <b>Reconstruction</b> plots ...</li></ul>
        """
        rec_time = round(reconstruction_result.time; digits=3)
        toast!(w, 2, "Reconstruction successful<br>Time: $rec_time s", body)
        w.img[] = reconstruction_result.value
    catch error
        @error "Reconstruction failed" exception=(error, catch_backtrace())
        restore_content!(w, previous_content, previous_state)
        failure_toast!(w, 2, "Reconstruction", error)
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

function update_bonito_progress!(w::KomaWindow, block, Nblocks)
    progress = floor(Int, block / Nblocks * 100)
    evaljs(w, js"""
        const bar = document.getElementById('simul_progress');
        if (bar) {
            bar.closest('button').style.setProperty('--simulation-progress', $progress + '%');
            bar.textContent = $progress + '%';
            bar.setAttribute('aria-valuenow', $progress);
        }
    """)
    return nothing
end

function ui_progressbar_callback(w::KomaWindow)
    return Callback(1, (progress_info, _, _, _) -> begin
        update_bonito_progress!(w, progress_info.block, progress_info.Nblocks)
    end)
end
