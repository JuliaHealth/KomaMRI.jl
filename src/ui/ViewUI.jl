const AcousticBand = @NamedTuple{freq::Float64, bw::Float64}

function read_acoustic_bands(filename)
    frequencies, bandwidths = Dict{Int,Float64}(), Dict{Int,Float64}()
    stem, extension = splitext(filename)
    companion = stem * "_GSWD_SAFETY" * extension
    for path in (isfile(companion) ? (filename, companion) : (filename,))
        for line in eachline(path)
            strip(line) == "### ASCCONV END ###" && break
            entry = match(r"^\s*asGPAParameters\[0\]\.sGCParameters\.aflAcousticResonance(Frequency|Bandwidth)\[(\d+)\]\s*=\s*([^\s/#]+)", line)
            isnothing(entry) && continue
            values = entry[1] == "Frequency" ? frequencies : bandwidths
            values[parse(Int, entry[2])] = parse(Float64, entry[3])
        end
    end
    isempty(frequencies) && throw(ArgumentError("No acoustic resonance frequencies found in ASC file."))
    return AcousticBand[(freq=frequencies[i], bw=bandwidths[i])
        for i in sort!(collect(keys(frequencies))) if frequencies[i] > 0]
end

function spectrum_gradient_samples(gradient, dt)
    samples = KomaMRIBase.event_samples(gradient)
    if isempty(samples.t) && dur(gradient) > gradient.delay
        return (t=[gradient.delay, dur(gradient)], A=[0.0, 0.0])
    end
    return (t=samples.t, A=samples.A .* γ)
end

function spectrum_gradient_samples(gradient::Union{KomaMRIBase.UniformlySampledGrad,KomaMRIBase.TimeShapedGrad}, dt)
    samples = KomaMRIBase.event_samples(gradient)
    isempty(samples.t) && return (t=[gradient.delay, dur(gradient)], A=[0.0, 0.0])
    amplitudes = gradient.A .* γ
    n = length(amplitudes)
    raster_centers = all(i -> abs((samples.t[i + 1] - gradient.delay) / dt - (i - 0.5)) < 1e-6, 1:n)
    raster_centers || return (t=samples.t, A=samples.A .* γ)

    # Match Pulseq's restoreAdditionalShapeSamples, including its Hz/m thresholds.
    restored = [gradient.first * γ; zeros(n)]
    for i in 1:n
        restored[i + 1] = 2amplitudes[i] - restored[i]
    end
    scale = maximum(abs, amplitudes)
    abs(last(restored) - gradient.last * γ) > 2e-5 * scale &&
        return (t=samples.t, A=samples.A .* γ)
    interpolated = [gradient.first * γ; (amplitudes[1:end-1] .+ amplitudes[2:end]) ./ 2; gradient.last * γ]
    edges = ifelse.(abs.(restored .- interpolated) .<= eps() + 2e-5 * scale, interpolated, restored)
    waveform = zeros(2n + 1)
    waveform[1:2:end] .= edges
    waveform[2:2:end] .= amplitudes
    changes = abs.([1; diff(diff(waveform)); 1]) .> 1e-8
    return (t=(gradient.delay .+ (0:2n) .* dt ./ 2)[changes], A=waveform[changes])
end

function spectrum_waveforms(seq, dt)
    starts = get_block_start_times(seq)
    return map(axes(seq.GR, 1)) do axis
        t, A = Float64[], Float64[]
        for block in eachindex(seq.DUR)
            gradient = seq.GR[axis, block]
            samples = spectrum_gradient_samples(gradient, dt)
            isempty(samples.t) && continue
            times = starts[block] .+ samples.t
            amplitudes = samples.A
            # Pulseq inserts half-raster ramps when nonzero gradients are separated by a gap.
            if !isempty(t) && last(t) + dt < first(times)
                if abs(last(A)) > 1e-6
                    push!(t, last(t) + dt / 2)
                    push!(A, 0)
                else
                    A[end] = 0
                end
                if abs(first(amplitudes)) > 1e-6
                    pushfirst!(times, first(times) - dt / 2)
                    pushfirst!(amplitudes, 0)
                else
                    amplitudes[1] = 0
                end
            end
            for (time, amplitude) in zip(times, amplitudes)
                if isempty(t) || time > last(t)
                    push!(t, time)
                    push!(A, amplitude)
                end
            end
        end
        (; t, A)
    end
end

function gradient_spectrum(seq; fmax=3000)
    dt = seq.DEF["GradientRasterTime"]
    nwin, oversampling = 5000, 3
    frequencies = (0:(nwin ÷ 2 - 1)) ./ nwin ./ dt ./ oversampling
    nf = count(<=(fmax), frequencies)
    waveforms = spectrum_waveforms(seq, dt)
    tmax = maximum(samples -> isempty(samples.t) ? 0.0 : last(samples.t), waveforms)
    tmax > 0 || throw(ArgumentError("Empty sequence passed to gradient spectrum."))
    nt = ceil(Int, tmax / dt)
    nsegments = cld(nt, nwin)
    sample_times = ((1:nt) .- 0.5) .* dt
    window = 0.5 .* (1 .- cos.(2π .* (1:nwin) ./ nwin))
    Rax = zeros(length(waveforms), nf)
    for (axis, samples) in enumerate(waveforms)
        waveform = KomaMRIBase.linear_interpolate_samples(samples, sample_times)
        append!(waveform, zeros(nsegments * nwin - nt))
        for start in 1:(nwin ÷ 2):(length(waveform) - nwin + 1)
            segment = zeros(nwin * oversampling)
            segment[1:nwin] .= @view waveform[start:(start + nwin - 1)]
            # Pulseq removes DC across the padded window, before applying the Hann window.
            segment .-= sum(segment) / length(segment)
            nsegments > 1 && (segment[1:nwin] .*= window)
            Rax[axis, :] .+= abs2.(@view fft(segment)[1:nf])
        end
    end
    Rax .= sqrt.(Rax ./ (2nsegments - 1))
    return (; R=vec(sqrt.(sum(abs2, Rax; dims=1))), Rax, F=frequencies[1:nf])
end

function plot_spectrum(seq, bands; darkmode=false)
    (; R, Rax, F) = gradient_spectrum(seq)
    data = [PlotlyBase.scatter(; x=F, y=Rax[axis, :], mode="lines", name)
        for (axis, name) in enumerate(("Gx", "Gy", "Gz"))]
    push!(data, PlotlyBase.scatter(; x=F, y=R, mode="lines", name="Gtot"))
    band_color = "#b57e83"
    band_regions = [PlotlyBase.attr(; type="rect", x0=band.freq - band.bw / 2,
        x1=band.freq + band.bw / 2, yref="paper", y0=0, y1=1,
        fillcolor=band_color, opacity=0.18, line_width=0, layer="below") for band in bands]
    band_lines = [PlotlyBase.attr(; type="line", x0=frequency, x1=frequency,
        yref="paper", y0=0, y1=1, line=PlotlyBase.attr(; color=band_color, dash))
        for band in bands for (frequency, dash) in
            ((band.freq, "solid"), (band.freq - band.bw / 2, "dash"), (band.freq + band.bw / 2, "dash"))]
    background, text, plot_background, grid, _ = KomaMRIPlots.theme_chooser(darkmode)
    layout = PlotlyBase.Layout(; xaxis=PlotlyBase.attr(; title="Frequency / Hz", range=[0, 3000], gridcolor=grid),
        yaxis=PlotlyBase.attr(; title="Gradient spectrum / Hz/m", rangemode="tozero", gridcolor=grid,
            exponentformat="power"),
        legend=PlotlyBase.attr(; orientation="h"), shapes=[band_regions; band_lines],
        annotations=[PlotlyBase.attr(; x=band.freq, y=1, yref="paper", yanchor="bottom",
            text="$(band.freq) Hz", font_color=band_color, showarrow=false) for band in bands],
        paper_bgcolor=background, plot_bgcolor=plot_background, font_color=text,
        margin=PlotlyBase.attr(; t=40, r=24, b=50, l=70))
    return PlotlyBase.Plot(data, layout; config=PlotlyBase.PlotConfig(; displaylogo=false, responsive=true))
end

function show_sequence!(w, seq, view; darkmode=w.darkmode[], physio=NoPhysioSignal())
    if view === :sequence
        display_loading!(w, "Plotting sequence ...")
        time_end = length(seq) > 1_000 ? dur(seq) * 1e3 : 30
        plot = plot_seq(
            seq;
            darkmode,
            range=[0 time_end],
            adaptive=true,
            slider=true,
            physio,
        )
        set_content!(w, plot_node(plot), "sequence")
    elseif view === :kspace
        display_loading!(w, "Plotting kspace ...")
        plot = plot_kspace(seq; darkmode, adaptive=true, adaptive_navigator=true, physio)
        set_content!(w, plot_node(plot), "kspace")
    elseif view in (:moment0, :moment1, :moment2)
        index = findfirst(==(view), (:moment0, :moment1, :moment2))
        order = index - 1
        display_loading!(w, "Plotting moment $order ...")
        plot = (plot_M0, plot_M1, plot_M2)[index](seq; adaptive=true, slider=true, darkmode, physio)
        set_content!(w, sequence_comparison(seq, plot; darkmode, physio, show_rf_center=true), "m$order")
    elseif view === :slew_rate
        display_loading!(w, "Plotting slew rate ...")
        plot = plot_slew_rate(seq; adaptive=true, slider=true, darkmode, physio)
        set_content!(w, sequence_comparison(seq, plot; darkmode, physio), "slew_rate")
    elseif view === :spectrum
        display_loading!(w, "Plotting gradient spectrum ...")
        plot = fetch(Threads.@spawn plot_spectrum(resolve_triggers(seq, physio), w.forbidden_bands[]; darkmode))
        set_content!(w, plot_node(plot), "spectrum")
    else
        throw(ArgumentError("Unsupported sequence view: $view"))
    end
    return nothing
end

function sequence_comparison(seq, plot; darkmode, physio, show_rf_center=false)
    waveform = plot_seq(seq; adaptive=true, slider=false, darkmode, physio)
    template = KomaMRIPlots.plot_template(waveform.source)
    if show_rf_center
        for trace in template.data
            get(trace, :name, "") == "RF_center" && (trace[:visible] = true)
        end
    end
    layout = template.layout
    layout[:updatemenus] = []
    layout[:xaxis][:showticklabels] = false
    layout[:margin][:b] = 6
    return DOM.div(waveform, plot;
        class="koma-sequence-comparison", dataLinkedTimePlots=true)
end

function show_phantom!(w, obj; key=:ρ, darkmode=w.darkmode[], slice_view=false)
    display_loading!(w, "Plotting phantom ...")
    plotter = slice_view ? plot_phantom_slice : plot_phantom
    plot = plotter(obj; key, time_samples=5, darkmode, adaptive=true, height=nothing)
    if !slice_view
        only(plot.template.layout[:sliders])[:pad][:b] = 4
    end
    return set_content!(w, plot_node(plot), slice_view ? "phantom_slice" : "phantom")
end

function show_scanner!(w, sys; darkmode=w.darkmode[])
    display_loading!(w, "Plotting receive sensitivities ...")
    return set_content!(w, plot_node(plot_coil_sens(sys; darkmode, adaptive=true, height=nothing)), "coils")
end

function show_scanner_parameters!(w, sys)
    display_loading!(w, "Displaying hardware limits ...")
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
    return set_content!(w, dictionary_page(values, "Hardware limits"), "scanneparams")
end

function observe_scanner!(w, scanner)
    receiver = scanner[].receiver
    limits = Tuple(getfield(scanner[].limits, name) for name in fieldnames(HardwareLimits))
    push!(w.listeners, on(scanner) do sys
        next_limits = Tuple(getfield(sys.limits, name) for name in fieldnames(HardwareLimits))
        receiver_changed = sys.receiver !== receiver
        limits_changed = !isequal(next_limits, limits)
        receiver, limits = sys.receiver, next_limits
        if receiver_changed
            show_scanner!(w, sys)
        elseif limits_changed || w.state[] == "scanneparams"
            show_scanner_parameters!(w, sys)
        else
            show_scanner!(w, sys)
        end
        return nothing
    end)
    return nothing
end

function show_parameters!(w, parameters, title, state)
    display_loading!(w, "Displaying $(lowercase(title)) ...")
    return set_content!(w, dictionary_page(parameters, title), state)
end

function show_configuration!(w)
    set_content!(w, DOM.div(
        dictionary_page(w.sim_params[], "Simulation [sim_params]"),
        dictionary_page(w.rec_params[], "Reconstruction [recon_params]");
        class="koma-configuration-view"), "configuration")
    return nothing
end

function show_signal!(w, raw; darkmode=w.darkmode[])
    display_loading!(w, "Plotting raw signal ...")
    return set_content!(w, plot_node(plot_signal(raw; adaptive=true, slider=true, darkmode)), "sig")
end

const KSPACE_DYNAMIC_RANGE_DB = 60

function image_values(img, type; limits=nothing)
    array = Array(img)
    data, zmin, zmax = if type === :absi
        values = array * prod(size(array)[1:2])
        zmin, zmax = isnothing(limits) ? extrema(abs, values) : limits
        values, zmin, zmax
    elseif type === :absk
        dynamic_range_db = KSPACE_DYNAMIC_RANGE_DB
        magnitude = abs.(fftc(array; dims=1:min(3, ndims(array))))
        peak = maximum(magnitude)
        floor = oftype(peak, 10^(-dynamic_range_db / 20))
        values = if iszero(peak)
            fill(oftype(peak, -dynamic_range_db), size(magnitude))
        else
            oftype(peak, 20) .* log10.(max.(magnitude ./ peak, floor))
        end
        values, -dynamic_range_db, 0
    else
        throw(ArgumentError("Unsupported image view: $type"))
    end
    return data, zmin, zmax
end

image_aspect_ratio(img) = 1
function image_aspect_ratio(img::AxisArrays.AxisArray)
    x, y = AxisArrays.axes(img)[1:2]
    return Float64(step(y.val) / step(x.val))
end

image_spacing(img) = (1,1,1)
image_spacing(img::AxisArrays.AxisArray) = ntuple(d -> step(AxisArrays.axes(img)[d].val),3)

function reconstruction_plot(data, type, darkmode; zmin, zmax, aspect_ratio=1, spacing=(1,1,1))
    plot = ndims(data)==3 ? plot_image(data; zmin,zmax,darkmode,spacing) : plot_image(data; zmin,zmax,darkmode)
    ndims(data)==3 && return plot
    plot.layout.margin[:t] = 6
    if size(data, 2) != 1
        plot.layout.yaxis[:constrain] = "domain"
        plot.layout.yaxis[:scaleratio] = aspect_ratio
        for trace in plot.data
            colorbar = get!(trace.fields, :colorbar, PlotlyBase.attr())
            colorbar[:ypad] = 0
        end
    elseif type === :absk
        plot.layout.yaxis[:title] = "Magnitude (dB)"
    elseif length(plot.data) == 1
        plot.layout.yaxis[:range] = [zmin, 1.1 * zmax]
    end
    return plot
end

function image_plot(img::AbstractArray, type, darkmode)
    data, zmin, zmax = image_values(img, type)
    aspect_ratio = type === :absi ? image_aspect_ratio(img) : 1
    values = size(img,3)==1 ? data[:,:,1] : data
    return plot_node(reconstruction_plot(values,type,darkmode; zmin,zmax,aspect_ratio,
        spacing=type===:absi ? image_spacing(img) : (1,1,1)); fit_colorbar=ndims(values)==2)
end

function reconstruction_slider(label, values)
    slider = Slider(eachindex(values[]); class="koma-reconstruction-slider-input")
    on(values) do options
        length(slider.values[]) == length(options) && return
        slider.values[] = collect(eachindex(options))
        slider.index[] = min(slider.index[], length(options))
    end
    selected = map((options, index) -> options[index], values, slider.index)
    tick_marks = map(values) do options
        DOM.span(
            (DOM.span(; class="koma-reconstruction-slider-tick") for _ in options)...;
            class="koma-reconstruction-slider-ticks",
        )
    end
    value_label = map(label, selected)
    style = map(options -> length(options) > 1 ? "" : "display: none", values)
    control = DOM.label(
        DOM.span(
            slider,
            tick_marks;
            class="koma-reconstruction-slider-track",
        ),
        DOM.span(value_label; class="koma-reconstruction-slider-label");
        class="koma-reconstruction-slider",
        style,
    )
    return selected, slider, control
end

function reconstruction_dropdown(name, values)
    dropdown = Dropdown(
        values;
        option_to_string=value -> uppercasefirst(replace(string(value), "_" => " ")),
        style=nothing,
        class="koma-reconstruction-dropdown-input",
    )
    control = DOM.label(
        dropdown,
        DOM.span(string(name); class="koma-reconstruction-slider-label");
        class="koma-reconstruction-dropdown",
    )
    selected = map(index -> values[index], dropdown.option_index)
    return selected, control
end

function matching_reconstructions(images, names, values)
    return filter(images) do image
        all(getproperty(image.labels, name) == value for (name, value) in zip(names, values))
    end
end

function reconstruction_label_options(images, name, selections)
    names = first.(selections)
    observables = last.(selections)
    isempty(observables) && return Observable(unique(getproperty(image.labels, name) for image in images))
    return map(observables...) do values...
        selected = matching_reconstructions(images, names, values)
        unique(getproperty(image.labels, name) for image in selected)
    end
end

const _RECONSTRUCTION_AXIS_LABELS = (
    z=:PAR,
    echos=:ECO,
    coils=:COIL,
    repetitions=:REP,
)

function reconstruction_dimensions(entry)
    names = AxisArrays.axisnames(entry.image)
    image_axes = AxisArrays.axes(entry.image)
    return [
        begin
            name = getproperty(_RECONSTRUCTION_AXIS_LABELS, names[dimension])
            values = if names[dimension] === :z
                # Image planes, not acquired k-space partition counters.
                collect(0:(size(entry.image, dimension) - 1))
            else
                source = get(entry.source, name, ())
                collect(length(source) == size(entry.image, dimension) ?
                    source : image_axes[dimension].val)
            end
            name => values
        end for dimension in 3:ndims(entry.image)
        if hasproperty(_RECONSTRUCTION_AXIS_LABELS, names[dimension])
    ]
end

function reconstruction_dimension(entry, name)
    return last(only(dimension for dimension in reconstruction_dimensions(entry) if first(dimension) === name))
end

function reconstruction_limits(selected, echo, type)
    type === :absk && return (-Float64(KSPACE_DYNAMIC_RANGE_DB), 0.0)
    return selected.magnitude_limits[findfirst(==(echo), selected.source.ECO)]
end

function reconstruction_series_plot(reconstruction, type, darkmode)
    reconstruction.policy.COIL === :rss && type !== :absi && return DOM.div(
        "RSS reconstruction contains magnitude only.";
        class="koma-reconstruction-message",
    )
    images = reconstruction.images
    stack_labels = filter(!=(:SLC), keys(first(images).labels))
    stack_slices = length(unique(image.labels.SLC for image in images))>1 && all(images) do image
        peers = matching_reconstructions(images, stack_labels, getproperty.(Ref(image.labels), stack_labels))
        size(image.image,3)==1 && all(peers) do peer
            size(peer.image)==size(image.image) && image_aspect_ratio(peer.image)==image_aspect_ratio(image.image) &&
                reconstruction_dimensions(peer)==reconstruction_dimensions(image)
        end
    end
    label_names = filter(keys(first(images).labels)) do name
        !(stack_slices && name===:SLC) && length(unique(getproperty(image.labels,name) for image in images))>1
    end
    label_selections = Pair{Symbol,Observable}[]
    controls = Any[]
    for name in label_names
        options = reconstruction_label_options(images, name, label_selections)
        selected, control = if name === :ROLE
            reconstruction_dropdown(name, options[])
        else
            selected, _, control = reconstruction_slider(options) do value
                "$(name): $value"
            end
            selected, control
        end
        push!(label_selections, name => selected)
        push!(controls, control)
    end

    selected_images = if isempty(label_selections)
        Observable(images)
    else
        selected_labels = first.(label_selections)
        map(last.(label_selections)...) do values...
            matching_reconstructions(images, selected_labels, values)
        end
    end
    selected_image = map(first,selected_images)
    dimension_names = filter(first.(reconstruction_dimensions(first(images)))) do name
        name!==:PAR && any(length(reconstruction_dimension(image,name))>1 for image in images)
    end
    dimension_selections = Pair{Symbol,Observable}[]
    dimension_values = Pair{Symbol,Observable}[]
    for name in dimension_names
        options = map(image -> reconstruction_dimension(image, name), selected_image)
        selected, slider, control = reconstruction_slider(options) do value
            "$(uppercase(string(name))): $value"
        end
        push!(dimension_selections, name => slider.index)
        push!(dimension_values, name => selected)
        push!(controls, control)
    end

    echo_index = findfirst(pair -> first(pair) === :ECO, dimension_values)
    echo = isnothing(echo_index) ? map(image -> only(image.source.ECO), selected_image) :
        last(dimension_values[echo_index])
    image = map(selected_images) do entries
        ordered = sort(entries;by=entry -> entry.labels.SLC)
        data = [type===:absk ? first(image_values(entry.image,type)) : entry.image for entry in ordered]
        (; entry=first(ordered),data=stack_slices ? cat(data...;dims=3) : only(data),
            zlabels=stack_slices ? [entry.labels.SLC for entry in ordered] : nothing)
    end
    limits = map(selected_image, echo) do image, echo
        reconstruction_limits(image, echo, type)
    end
    plot = map(image, limits, last.(dimension_selections)...) do image, limits, selected...
        selection = Dict(first.(dimension_selections) .=> selected)
        axis_names = AxisArrays.axisnames(image.entry.image)
        indices = ntuple(ndims(image.entry.image)) do dimension
            dimension <= 3 && return Colon()
            name = getproperty(_RECONSTRUCTION_AXIS_LABELS, axis_names[dimension])
            return get(selection, name, 1)
        end
        volume = @view image.data[indices...]
        plane = size(volume,3)==1 ? volume[:,:,1] : volume
        data = type === :absk ? plane : first(image_values(plane, type; limits))
        figure = reconstruction_plot(
            data, type, darkmode;
            zmin=first(limits),
            zmax=last(limits),
            aspect_ratio=type === :absi ? image_aspect_ratio(image.entry.image) : 1,
            spacing=type!==:absi ? (1,1,1) : stack_slices ?
                (1,image_aspect_ratio(image.entry.image),1) : image_spacing(image.entry.image),
        )
        if !isnothing(image.zlabels) && ndims(plane)==3
            model=figure.layout[:meta]["image_slices"]
            figure.layout[:meta]["image_slices"]=(;model...,zlabels=image.zlabels)
        end
        figure
    end
    return DOM.div(
        DOM.div(controls...; class="koma-reconstruction-controls"),
        plot_node(plot; fit_colorbar=true);
        class="koma-plot-stack",
    )
end

function image_plot(reconstruction::ReconstructionResult, type, darkmode)
    return reconstruction_series_plot(reconstruction, type, darkmode)
end

function show_image!(w, img, view; darkmode=w.darkmode[])
    messages = Dict(
        :absi => "Plotting image ...",
        :absk => "Plotting image k ...",
    )
    display_loading!(w, messages[view])
    return set_content!(w, image_plot(img, view, darkmode), string(view))
end
