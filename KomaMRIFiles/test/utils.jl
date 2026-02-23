"""
    quantize_to_pulseq(seq::Sequence)

Recreates exactly what `read_seq(write_seq(seq, ...))` would return: applies the same
quantization that write (register_adc!, etc.) and read (get_ADC, block DUR) apply, so that
`quantize_to_pulseq(seq) ≈ read_seq(write_seq(seq, filename))`.
"""
function quantize_to_pulseq(seq::Sequence)
    blockDurationRaster = KomaMRIFiles.get_blockDurationRaster(seq)
    qseq = deepcopy(seq)
    for bi in 1:length(qseq)
        # RFs: 1 µs raster (same as write)
        qseq.RF[1, bi].delay = quantize_time(qseq.RF[1, bi].delay, 1e-6)
        # Gradients: 1 µs raster
        for gi in 1:3
            qseq.GR[gi, bi].T     = quantize_time(qseq.GR[gi, bi].T, 1e-6)
            qseq.GR[gi, bi].delay = quantize_time(qseq.GR[gi, bi].delay, 1e-6)
            qseq.GR[gi, bi].rise  = quantize_time(qseq.GR[gi, bi].rise, 1e-6)
            qseq.GR[gi, bi].fall  = quantize_time(qseq.GR[gi, bi].fall, 1e-6)
        end
        # ADCs: mirror write (register_adc!) → file → read (get_ADC)
        # Write: dwell_ns = round(dwell_s*1e9), del_us = round((delay - dwell_s_exact/2)*1e6)
        # Read: dwell = dwell_ns*1e-9, delay = del_us*1e-6 + dwell/2, T = (N-1)*dwell
        adc = qseq.ADC[bi]
        dwell_s_raw = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
        dwell_s = round(Int, dwell_s_raw * 1e9) * 1e-9
        qseq.ADC[bi].T = adc.N == 1 ? dwell_s : (adc.N - 1) * dwell_s
        del_us = round(Int, (adc.delay - dwell_s/2) * 1e6)
        del_us = max(0, del_us)
        qseq.ADC[bi].delay = del_us * 1e-6 + dwell_s/2
        # Block DUR: write uses max(dur(block), dur(events)) then ceil to raster; 
        #read gives DUR = duration*raster. Use original seq's GR/RF durations (as in write) and simulated ADC end time.
        adc_end = qseq.ADC[bi].delay + qseq.ADC[bi].T - dwell_s/2
        max_event_duration = max(dur.(seq.GR[:, bi])..., dur(seq.RF[1, bi]), adc_end)
        block_duration = max(seq.DUR[bi], max_event_duration)
        qseq.DUR[bi] = ceil(block_duration / blockDurationRaster) * blockDurationRaster
    end
    return qseq
end

quantize_time(t::Array, factor) = t
quantize_time(t::Number,factor) = round(Int, t / factor) * factor

function round_trip_sequences()
    sequences = Sequence[]

    # 1. FID sequence
    seq = Sequence()
    seq += RF(1.0, 100e-6)
    seq += ADC(100, 100e-3, 50e-3)
    seq.DEF["Name"] = "fid"
    push!(sequences, seq)

    # 2. EPI sequence
    seq = PulseDesigner.EPI_example()
    seq.DEF["Name"] = "epi"
    push!(sequences, seq)

    # 3. Uniformly-shaped gradient
    t = 0:0.25:1
    A = 10*10e-6 * sqrt.(π*t) .* sin.(π*t)
    T = 10e-3
    delay, rise, fall = 1e-3, 5e-6 , 5e-6
    first, last = 1e-6, -5e-6
    gr = Grad(A, T, rise, fall, delay, first, last)
    seq = Sequence([gr])
    seq.DEF["Name"] = "gr-uniformly-shaped"
    push!(sequences, seq)

    return sequences
end