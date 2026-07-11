function launch_ui(;
    darkmode=true,
    frame=true,
    phantom_mode="2D",
    sim=Dict{String,Any}(),
    rec=Dict{Symbol,Any}(),
    sys=nothing,
    seq=nothing,
    obj=nothing,
    verbose=true,
    return_window=false,
    show_window=true,
    dev_tools=false,
)
    foreach(clear, (seq_ui, obj_ui, sys_ui, raw_ui, img_ui))
    versions = join([
        "KomaMRI.jl v$(pkgversion(KomaMRI))",
        "KomaMRIBase.jl v$(pkgversion(KomaMRIBase))",
        "KomaMRICore.jl v$(pkgversion(KomaMRICore))",
        "KomaMRIFiles.jl v$(pkgversion(KomaMRIFiles))",
        "KomaMRIPlots.jl v$(pkgversion(KomaMRIPlots))",
    ], "\n")
    w = setup_bonito_window(; darkmode, frame, dev_tools, versions)

    fieldnames_obj = [fieldnames(Phantom)[5:end-3]...]
    widgets_button_obj = [Button(string(field); style=nothing, class="btn btn-dark btn-sm m-1") for field in fieldnames_obj]

    sys_default = isnothing(sys) ? setup_scanner() : sys
    seq_default = isnothing(seq) ? setup_sequence(sys_default) : seq
    obj_default = isnothing(obj) ? setup_phantom(; phantom_mode) : obj
    sys_ui[] = sys_default
    seq_ui[] = seq_default
    obj_ui[] = obj_default
    raw_ui[] = setup_raw()
    img_ui[] = [0.0im 0.; 0. 0.]
    verbose && @info "Loaded default UI inputs" scanner="sys_ui[]" sequence="seq_ui[]" phantom="obj_ui[]" raw="raw_ui[]" image="img_ui[]"

    sim_params = merge(Dict{String,Any}(), sim)
    rec_params = merge(Dict{Symbol,Any}(:reco => "direct"), rec)
    seq_file = Ref("")

    if !(haskey(sim_params, "gpu") && sim_params["gpu"] == false)
        KomaMRICore.print_devices()
    end

    is_first_sim = true
    is_first_rec = true

    handle(w, "index") do _
        set_content!(w, w.home[], "index")
    end
    handle(w, "pulses_seq") do _
        show_sequence!(w, seq_ui[], :sequence; darkmode)
    end
    handle(w, "reload_seq") do _
        isempty(seq_file[]) || (seq_ui[] = callback_filepicker(seq_file[], w, seq_ui[]))
    end
    handle(w, "pulses_kspace") do _
        show_sequence!(w, seq_ui[], :kspace; darkmode)
    end
    handle(w, "pulses_M0") do _
        show_sequence!(w, seq_ui[], :moment0; darkmode)
    end
    handle(w, "pulses_M1") do _
        show_sequence!(w, seq_ui[], :moment1; darkmode)
    end
    handle(w, "pulses_M2") do _
        show_sequence!(w, seq_ui[], :moment2; darkmode)
    end
    handle(w, "phantom") do _
        show_phantom!(w, obj_ui[], widgets_button_obj; key=:ρ, darkmode)
    end
    handle(w, "scanner") do _
        show_scanner!(w, sys_ui[])
    end
    handle(w, "sim_params") do _
        show_parameters!(w, sim_params, "Simulation parameters", "simparams")
    end
    handle(w, "sig") do _
        show_signal!(w, raw_ui[]; darkmode)
    end
    handle(w, "rec_params") do _
        show_parameters!(w, rec_params, "Reconstruction parameters", "recparams")
    end
    handle(w, "reconstruction_absI") do _
        show_image!(w, img_ui[], :absi; darkmode)
    end
    handle(w, "reconstruction_angI") do _
        show_image!(w, img_ui[], :angi; darkmode)
    end
    handle(w, "reconstruction_absK") do _
        show_image!(w, img_ui[], :absk; darkmode)
    end

    for type in ("all", "sequence", "phantom", "scanner", "raw", "image")
        event = Dict(
            "all" => "matfolder",
            "sequence" => "matfolderseq",
            "phantom" => "matfolderpha",
            "scanner" => "matfoldersca",
            "raw" => "matfolderraw",
            "image" => "matfolderima",
        )[type]
        handle(w, event) do _
            save_ui!(w, seq_ui[], obj_ui[], sys_ui[], raw_ui[], img_ui[], rec_params; type)
        end
    end

    handle(w, "simulate") do _
        message = is_first_sim ? "Precompiling and running simulation functions ..." : "Running simulation ..."
        is_first_sim = false
        display_loading!(w, message)
        start_simulation_progress!(w)

        raw_aux = try
            raw = simulate(obj_ui[], seq_ui[], sys_ui[]; sim_params, callbacks=(ui_progressbar_callback(w),))
            rawfile = joinpath(tempdir(), "Koma_signal.mrd")
            @info "Exporting to ISMRMRD file: $rawfile"
            save(ISMRMRDFile(rawfile), raw)
            raw
        finally
            finish_simulation_progress!(w)
        end

        sim_time = raw_aux.params["userParameters"]["sim_time_sec"]
        body = """
            <ul class="list-unstyled mb-0">
                <li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View raw signal" aria-label="View raw signal" onclick="KomaUI.notify('sig')"><i class="bi bi-search"></i></button> Updating <b>Raw signal</b> plots ...</li>
                <li><button type="button" class="btn btn-primary btn-circle btn-circle-sm m-1" title="Reconstruct" aria-label="Reconstruct" onclick="KomaUI.notify('recon')"><i class="bi bi-caret-right-fill"></i></button> Ready to <b>reconstruct</b>?</li>
            </ul>
        """
        update_filename!(w, "rawname", "Koma_signal.mrd")
        toast!(w, 1, "Simulation successful<br>Time: $sim_time s", body)
        raw_ui[] = raw_aux
    end

    handle(w, "recon") do _
        message = is_first_rec ? "Precompiling and running reconstruction functions ..." : "Running reconstruction ..."
        is_first_rec = false
        display_loading!(w, message)
        spinner = "<div class=\"spinner-border spinner-border-sm text-light\" role=\"status\"></div>"
        evaljs(w, js"document.getElementById('recon!').innerHTML = $(spinner);")

        raw_aux = raw_ui[]
        raw_aux.profiles = raw_aux.profiles[getproperty.(getproperty.(raw_aux.profiles, :head), :flags) .!= 268435456]
        acq_data = AcquisitionData(raw_aux)
        acq_data.traj[1].circular = false
        acq_data.traj[1].nodes = acq_data.traj[1].nodes[1:2, :] ./ maximum(2 * abs.(acq_data.traj[1].nodes[:]))
        Nx, Ny = raw_aux.params["reconSize"][1:2]
        rec_params[:reconSize] = (Nx, Ny)
        rec_params[:densityWeighting] = true

        @info "Running reconstruction ..."
        rec_aux = @timed reconstruction(acq_data, rec_params)
        image = reshape(rec_aux.value.data, Nx, Ny, :)
        evaljs(w, js"document.getElementById('recon!').innerHTML = 'Reconstruct!';")
        body = """
            <ul class="list-unstyled mb-0"><li><button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View reconstruction" aria-label="View reconstruction" onclick="KomaUI.notify('reconstruction_absI')"><i class="bi bi-search"></i></button> Updating <b>Reconstruction</b> plots ...</li></ul>
        """
        toast!(w, 2, "Reconstruction successful<br>Time: $(rec_aux.time) s", body)
        img_ui[] = image
    end

    on(seq -> show_sequence!(w, seq, :sequence; darkmode), seq_ui)
    on(obj -> show_phantom!(w, obj, widgets_button_obj; key=:ρ, darkmode), obj_ui)
    for (widget, key) in zip(widgets_button_obj, fieldnames_obj)
        on(widget.value) do _
            show_phantom!(w, obj_ui[], widgets_button_obj; key, darkmode)
        end
    end
    on(sys -> show_scanner!(w, sys), sys_ui)
    on(raw -> show_signal!(w, raw; darkmode), raw_ui)
    on(img -> show_image!(w, img, :absi; darkmode), img_ui)

    setup_filepickers!(w; seq_file)

    @info "KomaMRI loaded successfully 🚀" KomaMRI=string(pkgversion(KomaMRI)) KomaMRIBase=string(pkgversion(KomaMRIBase)) KomaMRICore=string(pkgversion(KomaMRICore)) KomaMRIFiles=string(pkgversion(KomaMRIFiles)) KomaMRIPlots=string(pkgversion(KomaMRIPlots))
    show_window && show!(w)
    return return_window ? w : nothing
end

function setup_filepickers!(w::BonitoWindow; seq_file=Ref(""))
    setup_filepicker!(
        w,
        "#seqfilepicker",
        "#seqname",
        ".seq (Pulseq)",
        seq_ui;
        accept=".seq,.seqk",
        selected_file=seq_file,
    )
    setup_filepicker!(
        w,
        "#phafilepicker",
        "#phaname",
        ".phantom (Koma)/.h5 (JEMRIS)",
        obj_ui;
        accept=".phantom,.h5",
    )
    setup_filepicker!(
        w,
        "#sigfilepicker",
        "#rawname",
        ".h5/.mrd (ISMRMRD)",
        raw_ui;
        accept=".h5,.mrd",
    )
    return nothing
end

upload_bytes(data) = UInt8.(data)
upload_bytes(data::MsgPack.Extension) = data.data

function setup_filepicker!(
    w::BonitoWindow,
    selector::String,
    current::String,
    label::String,
    output;
    accept,
    selected_file=nothing,
)
    upload = Observable{Any}(nothing)
    on(upload) do file
        isnothing(file) && return nothing
        name = String(get(file, "name", get(file, :name, "upload")))
        data = get(file, "data", get(file, :data, UInt8[]))
        bytes = upload_bytes(data)
        directory = mktempdir()
        filename = joinpath(directory, basename(name))
        open(filename, "w") do io
            write(io, bytes)
        end
        isnothing(selected_file) || (selected_file[] = filename)
        output[] = callback_filepicker(filename, w, output[])
        return nothing
    end
    push!(w.on_render, session -> Bonito.on_document_load(session, js"""
            const mount = () => {
                const host = document.querySelector($(selector));
                if (!host) {
                    requestAnimationFrame(mount);
                    return;
                }
                const input = document.createElement('input');
                input.type = 'file';
                input.accept = $(accept);
                input.style.display = 'none';
                const button = document.createElement('label');
                button.className = 'koma-file-input';
                button.title = $(label);
                const icon = document.createElement('i');
                icon.className = 'bi bi-folder2-open';
                const caption = document.createElement('span');
                caption.className = 'koma-file-name';
                const setCaption = name => {
                    caption.textContent = name.length > $(MAX_UI_FILENAME_CHARS)
                        ? name.slice(0, $(MAX_UI_FILENAME_CHARS - 3)) + '...'
                        : name;
                    caption.title = name;
                    button.title = name;
                };
                setCaption(document.querySelector($(current))?.innerText.trim() || $(label));
                button.append(icon, caption, input);
                host.replaceChildren(button);
                input.addEventListener('change', async event => {
                    const file = event.target.files[0];
                    if (!file) return;
                    setCaption(file.name);
                    const data = new Uint8Array(await file.arrayBuffer());
                    $(upload).notify({name: file.name, data: data});
                });
            };
            mount();
        """))
    return nothing
end

function start_simulation_progress!(w::BonitoWindow)
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

function finish_simulation_progress!(w::BonitoWindow)
    evaljs(w, js"""
        const button = document.getElementById('simulate!');
        button.textContent = 'Simulate!';
        button.disabled = false;
    """)
    return nothing
end

function update_bonito_progress!(w::BonitoWindow, block, Nblocks)
    progress = floor(Int, block / Nblocks * 100)
    evaljs(w, js"""
        const bar = document.getElementById('simul_progress');
        if (bar) {
            bar.style.width = $progress + '%';
            bar.innerHTML = $progress + '%';
            bar.setAttribute('aria-valuenow', $progress);
        }
    """)
    return nothing
end

function ui_progressbar_callback(w::BonitoWindow)
    return Callback(1, (progress_info, _, _, _) -> update_bonito_progress!(w, progress_info.block, progress_info.Nblocks))
end
