include("ui/ExportMATFunctions.jl")
include("ui/ExportUIFunctions.jl")

"""
This launches the Koma's UI
"""
function KomaUI(; darkmode=true, frame=true, phantom_mode="2D", sim=Dict{String,Any}(), rec=Dict{Symbol,Any}(), dev_tools=false, blink_show=true)

    # For phantom sub-buttons
    fieldnames_obj = [fieldnames(Phantom)[5:end-3]...]
    widgets_button_obj = button.(string.(fieldnames_obj))

    # Setup the Blink window
    w, index = setup_blink_window(; darkmode, frame, dev_tools, blink_show)

    # Setup default simulation inputs (they have observables)
    sys = setup_scanner()
    seq = setup_sequence(sys)
    obj = setup_phantom(; phantom_mode)
    raw = setup_raw()
    img =  [0.0im 0.; 0. 0.]

    # Define observables
    obs_sys = Observable{Scanner}(sys)
    obs_seq = Observable{Sequence}(seq)
    obs_obj = Observable{Phantom}(obj)
    obs_raw = Observable{RawAcquisitionData}(raw)
    obs_img = Observable{Array}(img)

    # Define parameters (they are just internal variables)
    sim_params = merge(Dict{String,Any}(), sim)
    rec_params = merge(Dict{Symbol,Any}(:reco=>"direct"), rec)
    mat_folder = tempdir()

    # Boleans to indicate first time for precompilation
    is_first_sim = true
    is_first_rec = true

    # Print GPU information
    @info "Loading GPUs"
    KomaMRICore.print_gpus()

    # Handle "View" sidebar buttons
    handle(w, "index") do _
        content!(w, "div#content", index)
        @js_ w document.getElementById("content").dataset.content = "index"
    end
    handle(w, "pulses_seq") do _
        view_ui!(obs_seq[], w; type="sequence", darkmode)
    end
    handle(w, "pulses_kspace") do _
        view_ui!(obs_seq[], w; type="kspace", darkmode)
    end
    handle(w, "pulses_M0") do _
        view_ui!(obs_seq[], w; type="moment0", darkmode)
    end
    handle(w, "pulses_M1") do _
        view_ui!(obs_seq[], w; type="moment1", darkmode)
    end
    handle(w, "pulses_M2") do _
        view_ui!(obs_seq[], w; type="moment2", darkmode)
    end
    handle(w, "phantom") do _
        view_ui!(obs_obj[], w, obs_seq[], widgets_button_obj; key=:ρ, darkmode)
    end
    handle(w, "scanner") do _
        view_ui!(obs_sys[], w)
    end
    handle(w, "sim_params") do _
        view_ui!(sim_params, w)
    end
    handle(w, "sig") do _
        view_ui!(obs_raw[], w; darkmode)
    end
    handle(w, "rec_params") do _
        view_ui!(rec_params, w)
    end
    handle(w, "reconstruction_absI") do _
        view_ui!(obs_img[], w; type="absi", darkmode)
    end
    handle(w, "reconstruction_angI") do _
        view_ui!(obs_img[], w; type="angi", darkmode)
    end
    handle(w, "reconstruction_absK") do _
        view_ui!(obs_img[], w; type="absk", darkmode)
    end

    # Handle "Save" to mat sidebar buttons
    handle(w, "matfolder") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="all")
    end
    handle(w, "matfolderseq") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="sequence")
    end
    handle(w, "matfolderpha") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="phantom")
    end
    handle(w, "matfoldersca") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="scanner")
    end
    handle(w, "matfolderraw") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="raw")
    end
    handle(w, "matfolderima") do _
        save_ui!(w, obs_seq[], obs_obj[], obs_sys[], obs_raw[], obs_img[], rec_params, mat_folder; type="image")
    end

    # Handle "Simulation" sidebar button
    handle(w, "simulate") do _

        # Display loading and disable simulation button
        msg = "Running simulation ..."
        if is_first_sim
            msg = "Precompiling and running simulation functions ..."
            is_first_sim = false
        end
        display_loading!(w, msg)
        @js_ w document.getElementById("simulate!").setAttribute("disabled", true)

        # Define progress bar in the simulation button
        progressbar = """
            <div class="progress" style="background-color: #27292d;">
            <div id="simul_progress" class="progress-bar" role="progressbar" style="width: 0%; transition:none;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">0%</div>
            </div>
        """
        @js_ w (@var progressbar = $progressbar; document.getElementById("simulate!").innerHTML=progressbar)

        # Perform simulation
        raw_aux = simulate(obs_obj[], obs_seq[], obs_sys[]; sim_params, w)

        # After simulation, display the text on the simulation button (not the progress bar)
        @js_ w document.getElementById("simulate!").innerHTML="Simulate!"

        # Save the raw signal to a file in the temporal directory
        rawfile = tempdir()*"/Koma_signal.mrd"
        @info "Exporting to ISMRMRD file: $rawfile"
        KomaMRICore.save(ISMRMRDFile(rawfile), raw_aux)

        # Display message on UI
        sim_time = raw_aux.params["userParameters"]["sim_time_sec"]
        @js_ w (
            @var sim_time = $sim_time;
            @var name = $(obs_obj[].name);
            document.getElementById("rawname").innerHTML = "Koma_signal.mrd";
            Toasty(
                "1",
                """Simulation successfull<br>Time: <a id="sim_time"></a> s""",
                """
                    <ul>
                        <li>
                            <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('sig', 1)"><i class="fa fa-search"></i></button>
                            Updating <b>Raw signal</b> plots ...
                        </li>
                        <li>
                            <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('recon', 1)"><i class="bi bi-caret-right-fill"></i></button>
                            Ready to <b>reconstruct</b>?
                        </li>
                    </ul>
                """
            );
            document.getElementById("sim_time").innerHTML = sim_time;
        )

        # Set the dataset content and enable simulation and reconstruction buttons
        @js_ w document.getElementById("content").dataset.content = "simulation"
        @js_ w document.getElementById("simulate!").removeAttribute("disabled");
        @js_ w document.getElementById("recon!").removeAttribute("disabled");

        # Update the value of the raw signal observable
        # this calls the view_ui to display the raw signal
        obs_raw[] = raw_aux
    end

    # Handle "Reconstruction" sidebar button
    handle(w, "recon") do _

        # Display loading
        msg = "Running reconstruction ..."
        if is_first_rec
            msg = "Precompiling and running reconstruction functions ..."
            is_first_rec = false
        end
        display_loading!(w, msg)

        # Update loading icon for button
        buffericon = """
            <div class="spinner-border spinner-border-sm text-light" role="status"></div>
        """
        @js_ w (@var buffericon = $buffericon; document.getElementById("recon!").innerHTML = buffericon)

        # Get the value of the raw signal and prepare for reconstruction
        raw_aux = obs_raw[]
        raw_aux.profiles = raw_aux.profiles[getproperty.(getproperty.(raw_aux.profiles, :head), :flags) .!= 268435456] #Extra profile in JEMRIS simulations
        acq_data = AcquisitionData(raw_aux)
        acq_data.traj[1].circular = false #Removing circular window
        acq_data.traj[1].nodes = acq_data.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acq_data.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
        Nx, Ny = raw_aux.params["reconSize"][1:2]
        rec_params[:reconSize] = (Nx, Ny)
        rec_params[:densityWeighting] = true

        # Perform reconstruction
        @info "Running reconstruction ..."
        rec_aux = @timed reconstruction(acq_data, rec_params)
        image  = reshape(rec_aux.value.data, Nx, Ny, :)

        # After Recon go to Image
        recon_time = rec_aux.time
        @js_ w (document.getElementById("recon!").innerHTML = "Reconstruct!")
        @js_ w (
            @var recon_time = $recon_time;
            Toasty(
                "2",
                """Reconstruction successfull<br>Time: <a id="recon_time"></a> s""",
                """
                    <ul>
                        <li>
                            <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('reconstruction_absI', 1)"><i class="fa fa-search"></i></button>
                            Updating <b>Reconstruction</b> plots ...
                        </li>
                    </ul>
                """
            );
            document.getElementById("recon_time").innerHTML = recon_time;
        )
        @js_ w document.getElementById("content").dataset.content = "reconstruction"

        # Update the value of the image observable
        # this calls the view_ui to display the image
        obs_img[] = image
    end

    # Define functionality of sequence filepicker widget
    widget_filepicker_seq = filepicker(".seq (Pulseq)/.seqk (Koma)"; accept=".seq,.seqk")
    map!((filename) -> callback_filepicker(filename, w, obs_seq[]), obs_seq, widget_filepicker_seq)
    on((seq) -> view_ui!(seq, w; type="sequence", darkmode), obs_seq)

    # Define functionality of phantom filepicker widget and sub-buttons
    widget_filepicker_obj = filepicker(".phantom (Koma)/.h5 (JEMRIS)"; accept=".phantom,.h5")
    map!((filename) -> callback_filepicker(filename, w, obs_obj[]), obs_obj, widget_filepicker_obj)
    on((obj) -> view_ui!(obj, w, obs_seq[], widgets_button_obj; key=:ρ, darkmode), obs_obj)
    for (widget_button, key) in zip(widgets_button_obj, fieldnames_obj)
        on((cnt) -> view_ui!(cnt, w, obs_obj[], obs_seq[], widgets_button_obj; key, darkmode), widget_button)
    end

    # Define functionality of raw filepicker widget
    widget_filepicker_raw = filepicker(".h5/.mrd (ISMRMRD)"; accept=".h5,.mrd")
    map!((filename) -> callback_filepicker(filename, w, obs_raw[]), obs_raw, widget_filepicker_raw)
    on((raw) -> view_ui!(raw, w; darkmode), obs_raw)

    # Define functionality when image observable changes (after reconstruction)
    on((img) -> view_ui!(img, w; type="absi", darkmode), obs_img)

    # Add filepicker widgets to the UI
    content!(w, "#seqfilepicker", widget_filepicker_seq, async=false)
    content!(w, "#phafilepicker", widget_filepicker_obj, async=false)
    content!(w, "#sigfilepicker", widget_filepicker_raw, async=false)

    # Update Koma version on UI
    version = string(KomaMRICore.__VERSION__)
    content!(w, "#version", version, async=false)
    @info "Currently using KomaMRICore v$version"

    # Return the Blink Window or not
    (dev_tools) && return w
    return nothing
end

"""
Auxiliary function to updates KomaUI's simulation progress bar.
"""
function update_blink_window_progress!(w::Window, block, Nblocks)
    progress = string(floor(Int, block / Nblocks * 100))
    @js_ w (@var progress = $progress;
    document.getElementById("simul_progress").style.width = progress + "%";
    document.getElementById("simul_progress").innerHTML = progress + "%";
    document.getElementById("simul_progress").setAttribute("aria-valuenow", progress))
    return nothing
end
