function run_simulation!(w, sim_params; initial=false)
    previous_content = w.content[]
    previous_state = w.state[]
    message = initial ?
        "Precompiling and running simulation functions ..." : "Running simulation ..."
    threads = get(sim_params, "Nthreads", Threads.nthreads())
    simulation_device = Ref("CPU ($threads thread$(threads == 1 ? "" : "s"))")
    details = get(sim_params, "gpu", true) ? "Selecting GPU backend ..." : simulation_device[]
    display_loading!(w, message; details)
    start_simulation_progress!(w)

    raw = try
        raw = simulate(
            obj_ui[],
            seq_ui[],
            sys_ui[];
            sim_params,
            callbacks=(ui_progressbar_callback(w, simulation_device),),
            physio=physio_ui[],
        )
        rawfile = joinpath(tempdir(), "Koma_signal.mrd")
        @info "Exporting to ISMRMRD file: $rawfile"
        save(ISMRMRDFile(rawfile), raw)
        raw
    catch error
        @error "Simulation failed" exception=(error, catch_backtrace())
        restore_content!(w, previous_content, previous_state)
        failure_toast!(w, 1, "Simulation", error)
        return nothing
    finally
        finish_simulation_progress!(w)
    end

    sim_time = round(raw.params["userParameters"]["sim_time_sec"]; digits=3)
    body = """
        <ul class="list-unstyled mb-0">
            <li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View raw signal" aria-label="View raw signal" onclick="KomaUI.notify('sig')"><i class="bi bi-search"></i></button> Updating <b>Raw signal</b> plots ...</li>
            <li><button type="button" class="btn btn-primary btn-circle btn-circle-sm m-1" title="Reconstruct" aria-label="Reconstruct" onclick="KomaUI.notify('recon')"><i class="bi bi-caret-right-fill"></i></button> Ready to <b>reconstruct</b>?</li>
        </ul>
    """
    update_filename!(w, "rawname", "Koma_signal.mrd")
    toast!(w, 1, "$(simulation_device[]) simulation successful<br>Time: $sim_time s", body)
    raw_ui[] = raw
    return nothing
end

function run_reconstruction!(w, rec_params; initial=false)
    previous_content = w.content[]
    previous_state = w.state[]
    message = initial ?
        "Precompiling and running reconstruction functions ..." : "Running reconstruction ..."
    display_loading!(w, message)
    spinner = "<div class=\"spinner-border spinner-border-sm text-light\" role=\"status\"></div>"
    evaljs(w, js"document.getElementById('recon!').innerHTML = $(spinner);")

    try
        raw = _imaging_raw_data(raw_ui[])
        acq_data = AcquisitionData(raw)
        acq_data.traj[1].circular = false
        acq_data.traj[1].nodes =
            acq_data.traj[1].nodes[1:2, :] ./ maximum(2 * abs.(acq_data.traj[1].nodes[:]))
        Nx, Ny = raw.params["reconSize"][1:2]
        rec_params[:reconSize] = (Nx, Ny)
        rec_params[:densityWeighting] = true

        @info "Running reconstruction ..."
        reconstruction_result = @timed reconstruction(acq_data, rec_params)
        image = reshape(reconstruction_result.value.data, Nx, Ny, :)
        body = """
            <ul class="list-unstyled mb-0"><li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View reconstruction" aria-label="View reconstruction" onclick="KomaUI.notify('reconstruction_absI')"><i class="bi bi-search"></i></button> Updating <b>Reconstruction</b> plots ...</li></ul>
        """
        rec_time = round(reconstruction_result.time; digits=3)
        toast!(w, 2, "Reconstruction successful<br>Time: $rec_time s", body)
        img_ui[] = image
    catch error
        @error "Reconstruction failed" exception=(error, catch_backtrace())
        restore_content!(w, previous_content, previous_state)
        failure_toast!(w, 2, "Reconstruction", error)
    finally
        evaljs(w, js"document.getElementById('recon!').innerHTML = 'Reconstruct!';")
    end
    return nothing
end

function start_simulation_progress!(w::KomaWindow)
    isnothing(w.session[]) && return nothing
    Bonito.evaljs_value(w.session[], js"""(() => {
        const button = document.getElementById('simulate!');
        button.disabled = true;
        const progress = document.createElement('div');
        progress.className = 'progress w-100';
        progress.style.backgroundColor = '#27292d';
        const bar = document.createElement('div');
        bar.id = 'simul_progress';
        bar.className = 'progress-bar';
        bar.style.width = '0%';
        bar.style.transition = 'none';
        bar.setAttribute('role', 'progressbar');
        bar.setAttribute('aria-valuenow', '0');
        bar.setAttribute('aria-valuemin', '0');
        bar.setAttribute('aria-valuemax', '100');
        bar.textContent = '0%';
        progress.appendChild(bar);
        button.replaceChildren(progress);
        return true;
    })()""")
    return nothing
end

function finish_simulation_progress!(w::KomaWindow)
    evaljs(w, js"""
        const button = document.getElementById('simulate!');
        button.textContent = 'Simulate!';
        button.disabled = false;
    """)
    return nothing
end

function update_bonito_progress!(w::KomaWindow, block, Nblocks, status)
    progress = floor(Int, block / Nblocks * 100)
    evaljs(w, js"""
        const bar = document.getElementById('simul_progress');
        if (bar) {
            bar.style.width = $progress + '%';
            bar.innerHTML = $progress + '%';
            bar.setAttribute('aria-valuenow', $progress);
        }
        const details = document.getElementById('loadstatus');
        if (details) details.textContent = $(status);
    """)
    return nothing
end

function ui_progressbar_callback(w::KomaWindow, simulation_device)
    gpu_backend = Ref{Union{Nothing,String}}(nothing)
    callback = function(progress_info, _, _, sim_params)
        threads = sim_params["Nthreads"]
        simulation_device[] = if sim_params["gpu"]
            if isnothing(gpu_backend[])
                gpu_backend[] = KomaMRICore.name(KomaMRICore.get_backend(true))
            end
            "GPU ($(gpu_backend[]))"
        else
            "CPU ($threads thread$(threads == 1 ? "" : "s"))"
        end
        method = nameof(typeof(sim_params["sim_method"]))
        status = "$(simulation_device[]) · $method · $(uppercase(sim_params["precision"]))"
        update_bonito_progress!(w, progress_info.block, progress_info.Nblocks, status)
    end
    return Callback(1, callback)
end
