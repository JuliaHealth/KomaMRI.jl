sequence_slider_visible(::NoPhysioSignal, long_seq) = !long_seq
sequence_slider_visible(::AbstractPhysioSignal, _) = false

function show_sequence!(w, seq, view; darkmode=true, physio=NoPhysioSignal())
    if view === :sequence
        display_loading!(w, "Plotting sequence ...")
        long_seq = length(seq) > 1_000
        time_end = long_seq ? dur(seq) * 1e3 : 30
        plot = plot_seq(
            seq;
            darkmode,
            range=[0 time_end],
            slider=sequence_slider_visible(physio, long_seq),
            gl=long_seq,
            show_adc=false,
            physio,
        )
        set_content!(w, plot_node(plot), "sequence")
    elseif view === :kspace
        display_loading!(w, "Plotting kspace ...")
        set_content!(w, plot_node(plot_kspace(seq; darkmode)), "kspace")
    elseif view === :moment0
        display_loading!(w, "Plotting moment 0 ...")
        set_content!(w, plot_node(plot_M0(seq; darkmode)), "m0")
    elseif view === :moment1
        display_loading!(w, "Plotting moment 1 ...")
        set_content!(w, plot_node(plot_M1(seq; darkmode)), "m1")
    elseif view === :moment2
        display_loading!(w, "Plotting moment 2 ...")
        set_content!(w, plot_node(plot_M2(seq; darkmode)), "m2")
    else
        throw(ArgumentError("Unsupported sequence view: $view"))
    end
    return nothing
end

function show_phantom!(w, obj, buttons; key=:ρ, darkmode=true)
    display_loading!(w, "Plotting phantom ...")
    plot = plot_phantom_map(obj, key; time_samples=5, darkmode)
    return set_content!(
        w,
        DOM.div(DOM.div(buttons...), plot_node(plot); class="koma-plot-stack"),
        "phantom",
    )
end

function show_scanner!(w, sys)
    display_loading!(w, "Displaying scanner parameters ...")
    values = [
        "B0" => sys.limits.B0,
        "B1" => sys.limits.B1,
        "Gmax" => sys.limits.Gmax,
        "Smax" => sys.limits.Smax,
        "ADC_dt" => sys.limits.ADC_Δt,
        "DUR_dt" => sys.limits.DUR_Δt,
        "GR_dt" => sys.limits.GR_Δt,
        "RF_dt" => sys.limits.RF_Δt,
        "RF_ring_down_time" => sys.limits.RF_ring_down_time,
        "RF_dead_time" => sys.limits.RF_dead_time,
        "ADC_dead_time" => sys.limits.ADC_dead_time,
    ]
    return set_content!(w, dictionary_page(values, "Scanner parameters"), "scanneparams")
end

function show_parameters!(w, parameters, title, state)
    display_loading!(w, "Displaying $(lowercase(title)) ...")
    return set_content!(w, dictionary_page(parameters, title), state)
end

function show_signal!(w, raw; darkmode=true)
    display_loading!(w, "Plotting raw signal ...")
    return set_content!(w, plot_node(plot_signal(raw; darkmode)), "sig")
end

function image_plot(img, type, darkmode)
    data, zmin, zmax = if type === :absi
        values = abs.(img) * prod(size(img)[1:2])
        values, minimum(values), maximum(values)
    elseif type === :angi
        values = angle.(img)
        values, -π, π
    elseif type === :absk
        values = log.(abs.(fftc(img)) .+ 1)
        values, 0, 0.1 * maximum(values)
    else
        throw(ArgumentError("Unsupported image view: $type"))
    end
    slider = Slider(1:size(img, 3))
    plot = map(slider.value) do slice
        title = "Reconstruction ($slice/$(size(img, 3)))"
        plot_node(plot_image(data[:, :, slice]; zmin, zmax, darkmode, title))
    end
    return size(img, 3) == 1 ? plot : DOM.div(slider, plot; class="koma-plot-stack")
end

function show_image!(w, img, view; darkmode=true)
    messages = Dict(
        :absi => "Plotting image magnitude ...",
        :angi => "Plotting image phase ...",
        :absk => "Plotting image k ...",
    )
    display_loading!(w, messages[view])
    return set_content!(w, image_plot(img, view, darkmode), string(view))
end
