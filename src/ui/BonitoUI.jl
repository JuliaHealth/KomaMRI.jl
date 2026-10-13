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
    browser=false,
    host="127.0.0.1",
    port=9384,
)
    versions = join([
        "KomaMRI.jl v$(pkgversion(KomaMRI))",
        "KomaMRIBase.jl v$(pkgversion(KomaMRIBase))",
        "KomaMRICore.jl v$(pkgversion(KomaMRICore))",
        "KomaMRIFiles.jl v$(pkgversion(KomaMRIFiles))",
        "KomaMRIPlots.jl v$(pkgversion(KomaMRIPlots))",
    ], "\n")
    sys_default = isnothing(sys) ? setup_scanner() : sys
    seq_default = isnothing(seq) ? setup_sequence(sys_default) : seq
    obj_default = isnothing(obj) ? setup_phantom(; phantom_mode) : obj
    sim_params = merge(Dict{String,Any}(), sim)
    rec_params = merge(Dict{Symbol,Any}(:reco => "direct", :recon_policy => DEFAULT_RECON_POLICY), rec)
    w = setup_bonito_window(
        sys_default, seq_default, obj_default, sim_params, rec_params;
        darkmode, frame, dev_tools, versions,
    )
    setup_filepickers!(w)
    setup_host_status!(w)
    setup_exports!(w)

    verbose && @info "Loaded default UI inputs" scanner="w.sys[]" sequence="w.seq[]" phantom="w.obj[]" physio="w.physio[]" raw="w.raw[]" image="w.img[]"

    if !(haskey(sim_params, "gpu") && sim_params["gpu"] == false)
        KomaMRICore.print_devices()
    end

    setup_actions!(w)
    push!(w.listeners, on(physio -> show_sequence!(w, w.seq[], :sequence; physio), w.physio))
    push!(w.listeners, on(seq -> w.physio[] = default_physio_signal(seq), w.seq))
    push!(w.listeners, on(obj -> show_phantom!(w, obj; slice_view=w.state[] == "phantom_slice"), w.obj))
    observe_scanner!(w, w.sys)
    push!(w.listeners, on(raw -> show_signal!(w, raw), w.raw))
    push!(w.listeners, on(img -> show_image!(w, img, :absi), w.img))
    push!(w.listeners, on(params -> w.state[] == "configuration" ? show_configuration!(w) :
        w.state[] == "simparams" && show_parameters!(w, params, "Simulation parameters", "simparams"), w.sim_params))
    push!(w.listeners, on(params -> w.state[] == "configuration" ? show_configuration!(w) :
        show_parameters!(w, params, "Reconstruction parameters", "recparams"), w.rec_params))

    if show_window
        if browser
            w.server[] = Bonito.Server(w.app, host, port)
            @info "Open KomaUI in your browser" url=Bonito.HTTPServer.online_url(w.server[], "/")
        else
            show!(w)
        end
    end
    @info "KomaMRI loaded successfully 🚀" KomaMRI=string(pkgversion(KomaMRI)) KomaMRIBase=string(pkgversion(KomaMRIBase)) KomaMRICore=string(pkgversion(KomaMRICore)) KomaMRIFiles=string(pkgversion(KomaMRIFiles)) KomaMRIPlots=string(pkgversion(KomaMRIPlots))
    return return_window ? w : nothing
end
