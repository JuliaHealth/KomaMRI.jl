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
    versions = join([
        "KomaMRI.jl v$(pkgversion(KomaMRI))",
        "KomaMRIBase.jl v$(pkgversion(KomaMRIBase))",
        "KomaMRICore.jl v$(pkgversion(KomaMRICore))",
        "KomaMRIFiles.jl v$(pkgversion(KomaMRIFiles))",
        "KomaMRIPlots.jl v$(pkgversion(KomaMRIPlots))",
    ], "\n")
    w = setup_bonito_window(; darkmode, frame, dev_tools, versions)
    seq_file = Ref("")
    setup_filepickers!(w; seq_file)
    show_window && show!(w)

    fieldnames_obj = [fieldnames(Phantom)[5:end-3]...]
    widgets_button_obj = [Button(string(field); style=nothing, class="btn btn-dark btn-sm m-1") for field in fieldnames_obj]

    sys_default = isnothing(sys) ? setup_scanner() : sys
    seq_default = isnothing(seq) ? setup_sequence(sys_default) : seq
    obj_default = isnothing(obj) ? setup_phantom(; phantom_mode) : obj
    physio_default = default_physio_signal(seq_default)
    sys_ui[] = sys_default
    seq_ui[] = seq_default
    obj_ui[] = obj_default
    physio_ui[] = physio_default
    raw_ui[] = setup_raw()
    img_ui[] = [0.0im 0.; 0. 0.]
    verbose && @info "Loaded default UI inputs" scanner="sys_ui[]" sequence="seq_ui[]" phantom="obj_ui[]" physio="physio_ui[]" raw="raw_ui[]" image="img_ui[]"

    sim_params = merge(Dict{String,Any}(), sim)
    rec_params = merge(Dict{Symbol,Any}(:reco => "direct"), rec)

    if !(haskey(sim_params, "gpu") && sim_params["gpu"] == false)
        KomaMRICore.print_devices()
    end

    is_first_sim = true
    is_first_rec = true

    handle(w, "index") do _
        set_content!(w, w.home[], "index")
    end
    handle(w, "pulses_seq") do _
        show_sequence!(w, seq_ui[], :sequence; darkmode, physio=physio_ui[])
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
        initial = is_first_sim
        is_first_sim = false
        run_simulation!(w, sim_params; initial)
    end
    handle(w, "recon") do _
        initial = is_first_rec
        is_first_rec = false
        run_reconstruction!(w, rec_params; initial)
    end
    push!(w.listeners, on(physio -> show_sequence!(w, seq_ui[], :sequence; darkmode, physio), physio_ui))
    push!(w.listeners, on(seq -> physio_ui[] = default_physio_signal(seq), seq_ui))
    push!(w.listeners, on(obj -> show_phantom!(w, obj, widgets_button_obj; key=:ρ, darkmode), obj_ui))
    for (widget, key) in zip(widgets_button_obj, fieldnames_obj)
        push!(w.listeners, on(widget.value) do _
            show_phantom!(w, obj_ui[], widgets_button_obj; key, darkmode)
        end)
    end
    push!(w.listeners, on(sys -> show_scanner!(w, sys), sys_ui))
    push!(w.listeners, on(raw -> show_signal!(w, raw; darkmode), raw_ui))
    push!(w.listeners, on(img -> show_image!(w, img, :absi; darkmode), img_ui))

    @info "KomaMRI loaded successfully 🚀" KomaMRI=string(pkgversion(KomaMRI)) KomaMRIBase=string(pkgversion(KomaMRIBase)) KomaMRICore=string(pkgversion(KomaMRICore)) KomaMRIFiles=string(pkgversion(KomaMRIFiles)) KomaMRIPlots=string(pkgversion(KomaMRIPlots))
    return return_window ? w : nothing
end
