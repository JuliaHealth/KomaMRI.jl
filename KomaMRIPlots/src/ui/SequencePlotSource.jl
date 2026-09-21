struct SequencePlotSource{B,T}
    blocks::B
    triggers::T
end

plot_template(source::SequencePlotSource) = source.blocks

const WAVEFORM_CHANNELS = (:Gx, :Gy, :Gz, :rf_abs, :rf_phase, :rf_frequency, :rf_center, :ADC)
const ADC_MARKER_SIZE = 6
const ADC_MARKER_WIDTH = 1

function block_samples(context, channel, block)
    (; seq, starts, labels, real_rf, ecg) = context
    text = Float64[]
    if channel === :ECG
        first_sample = searchsortedfirst(ecg[:x], starts[block] * 1e3)
        last_sample = block == length(seq) ? length(ecg[:x]) :
            searchsortedfirst(ecg[:x], starts[block + 1] * 1e3) - 1
        return (; x=ecg[:x][first_sample:last_sample], y=ecg[:y][first_sample:last_sample])
    elseif channel in (:Gx, :Gy, :Gz)
        samples = KomaMRIBase.event_samples(seq.GR[findfirst(==(channel), (:Gx, :Gy, :Gz)), block])
        t, y = samples.t, samples.A .* 1e3
    elseif channel in (:rf_abs, :rf_phase, :rf_frequency)
        rf = seq.RF[1, block]
        samples = channel === :rf_frequency ? KomaMRIBase.event_samples(rf, Val(:Δf)) : KomaMRIBase.event_samples(rf)
        t, amplitude = samples.t, samples.A
        y = channel === :rf_frequency ? amplitude .* 1e-3 :
            channel === :rf_abs ? (real_rf ? real.(amplitude) : abs.(amplitude)) .* 1e6 :
            (real_rf ? zeros(length(amplitude)) : angle.(amplitude))
        channel !== :rf_abs && (text = ones(length(t)))
    elseif channel === :rf_center
        rf = seq.RF[1, block]
        if is_RF_on(rf)
            t = [rf.delay + rf.center]
            value = first(first(KomaMRIBase.get_rfs(seq[block], t)))
            y, text = [abs(value) * 1e6], [angle(value)]
        else
            t = y = Float64[]
        end
    elseif channel === :ADC_samples
        samples = KomaMRIBase.event_samples(seq.ADC[block])
        t, y = samples.t, Float64.(samples.A)
    elseif channel === :ADC
        adc = seq.ADC[block]
        t = is_ADC_on(adc) ? (adc.N == 1 ? fill(adc.delay + adc.T / 2, 4) :
            [adc.delay, adc.delay, adc.delay + adc.T, adc.delay + adc.T]) : Float64[]
        y = isempty(t) ? Float64[] : [0.0, 1.0, 1.0, 0.0]
    else
        t = is_ADC_on(seq.ADC[block]) ? [seq.DUR[block] / 2] : Float64[]
        y = fill(Float64(getproperty(labels[block], channel)), length(t))
    end
    return (; x=(starts[block] .+ t) .* 1e3, y, text)
end

function sequence_source(seq; physio=NoPhysioSignal(), width=nothing, height=nothing,
    slider=true, darkmode=false, range=[], title="")
    seq = resolve_triggers(seq, physio)
    starts = get_block_start_times(seq)
    labels = get_labels(seq)
    real_rf = all(seq.RF) do rf
        samples = KomaMRIBase.event_samples(rf)
        all(amplitude -> isapprox(imag(amplitude), 0; atol=eps()), samples.A)
    end
    adc_blocks = findall(is_ADC_on, seq.ADC)
    label_names = [name for name in fieldnames(AdcLabels)
        if any(b -> !iszero(getproperty(labels[b], name)), adc_blocks)]
    channels = [WAVEFORM_CHANNELS...]
    template = plot_seq_nonadaptive(seq[1:1]; slider=false, show_adc=false, darkmode)
    data = template.data[1:length(WAVEFORM_CHANNELS)]
    samples = deepcopy(data[end])
    samples[:mode] = "markers"
    samples[:showlegend] = false
    samples[:marker][:size] = ADC_MARKER_SIZE
    samples[:marker][:symbol] = "line-ns"
    samples[:marker][:line] = Dict(:color=>samples[:marker][:color], :width=>ADC_MARKER_WIDTH)
    push!(data, samples)
    push!(channels, :ADC_samples)
    non_label_count = length(channels)
    append!(channels, label_names)
    for name in label_names
        trace = deepcopy(data[end])
        merge!(trace.fields, Dict(:name=>string(name), :legendgroup=>string(name),
            :showlegend=>false, :mode=>"markers", :visible=>false,
            :marker=>attr(color="white", symbol="x")))
        push!(data, trace)
    end
    layout, config = generate_seq_time_layout_config(
        title, width, height, range, slider, false, darkmode;
        T0=starts, label_to_show=label_names, non_label_count)
    physio_traces = typeof(first(data))[]
    _add_physio!(physio_traces, layout, scatter, seq, physio, "x")
    ecg = isempty(physio_traces) ? nothing : only(physio_traces)
    if !isnothing(ecg)
        push!(data, ecg)
        push!(channels, :ECG)
    end
    context = (; seq, starts, labels, real_rf, ecg)
    full_range = [0.0, dur(seq) * 1e3]
    break_blocks = [get(trace.fields, :mode, "lines") != "markers" for trace in data]
    !isnothing(ecg) && (break_blocks[end] = false)
    # Keep the original ECG arrays in the context; summarizing replaces plotted trace arrays.
    !isnothing(ecg) && (data[end] = deepcopy(ecg))
    source = summarize_blocks(context, channels, data, layout, config, starts .* 1e3, full_range, block_samples; break_blocks)
    trigger_times = [(starts[b] + delay) * 1e3 for b in eachindex(seq.DUR)
        for extension in seq.EXT[b] for delay in _trigger_delays(extension)]
    trigger_trace = trigger_line = nothing
    if !isempty(trigger_times)
        block = findfirst(b -> has_trigger(seq[b]), eachindex(seq.DUR))
        trigger_plot = plot_seq_nonadaptive(seq[block:block]; darkmode)
        trigger_trace = only(filter(trace -> get(trace.fields, :name, "") == "Trigger", trigger_plot.data))
        trigger_line = first(trigger_plot.layout[:shapes])
        layout[:yaxis3] = trigger_plot.layout[:yaxis3]
    end
    return SequencePlotSource(source, (; times=trigger_times, trace=trigger_trace, line=trigger_line))
end

function window_data(source::SequencePlotSource, range, width, visibility)
    (; seq, starts) = source.blocks.context
    first_block = clamp(searchsortedlast(starts, range[1] / 1e3), 1, length(seq))
    last_block = clamp(searchsortedlast(starts, range[2] / 1e3), 1, length(seq))
    min_spacing = ADC_MARKER_WIDTH * (range[2] - range[1]) / max(width - 70, 1)
    resolved, previous = true, -Inf
    for block in first_block:last_block
        times = (starts[block] .+ KomaMRIBase.times(seq.ADC[block])) .* 1e3
        inside = searchsortedfirst(times, range[1]):searchsortedlast(times, range[2])
        for time in @view(times[inside])
            if time - previous < min_spacing
                resolved = false
                break
            end
            previous = time
        end
        resolved || break
    end
    visibility = copy(visibility)
    visibility[length(WAVEFORM_CHANNELS) + 1] = visibility[length(WAVEFORM_CHANNELS)] === true && resolved
    payload = block_window_data(source.blocks, range, width, visibility)
    isempty(source.triggers.times) && return payload
    times = filter(t -> range[1] <= t <= range[2], source.triggers.times)
    layout = deepcopy(payload.layout)
    layout[:shapes] = [merge(copy(source.triggers.line.fields), Dict(:x0=>t, :x1=>t)) for t in times]
    trace = copy(source.triggers.trace.fields)
    positions = source.triggers.trace[:y]
    trace[:x] = repeat(times; inner=length(positions))
    trace[:y] = repeat(positions, length(times))
    trace[:uid] = "triggers"
    return (; payload..., data=[payload.data; [trace]], layout)
end
