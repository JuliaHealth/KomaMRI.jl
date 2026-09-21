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
    sys_default = isnothing(sys) ? setup_scanner() : sys
    seq_default = isnothing(seq) ? setup_sequence(sys_default) : seq
    obj_default = isnothing(obj) ? setup_phantom(; phantom_mode) : obj
    sim_params = merge(Dict{String,Any}(), sim)
    rec_params = merge(Dict{Symbol,Any}(:reco => "direct"), rec)
    w = setup_bonito_window(
        sys_default, seq_default, obj_default, sim_params, rec_params;
        darkmode, frame, dev_tools, versions,
    )
    setup_filepickers!(w)

    verbose && @info "Loaded default UI inputs" scanner="w.sys[]" sequence="w.seq[]" phantom="w.obj[]" physio="w.physio[]" raw="w.raw[]" image="w.img[]"

    if !(haskey(sim_params, "gpu") && sim_params["gpu"] == false)
        KomaMRICore.print_devices()
    end

    setup_actions!(w; darkmode)
    push!(w.listeners, on(physio -> show_sequence!(w, w.seq[], :sequence; darkmode, physio), w.physio))
    push!(w.listeners, on(seq -> w.physio[] = default_physio_signal(seq), w.seq))
    push!(w.listeners, on(obj -> show_phantom!(w, obj; darkmode), w.obj))
    observe_scanner!(w, w.sys; darkmode)
    push!(w.listeners, on(raw -> show_signal!(w, raw; darkmode), w.raw))
    push!(w.listeners, on(img -> show_image!(w, img, :absi; darkmode), w.img))
    push!(w.listeners, on(params -> show_parameters!(w, params, "Simulation parameters", "simparams"), w.sim_params))
    push!(w.listeners, on(params -> show_parameters!(w, params, "Reconstruction parameters", "recparams"), w.rec_params))

    show_window && show!(w)
    @info "KomaMRI loaded successfully 🚀" KomaMRI=string(pkgversion(KomaMRI)) KomaMRIBase=string(pkgversion(KomaMRIBase)) KomaMRICore=string(pkgversion(KomaMRICore)) KomaMRIFiles=string(pkgversion(KomaMRIFiles)) KomaMRIPlots=string(pkgversion(KomaMRIPlots))
    return return_window ? w : nothing
end
