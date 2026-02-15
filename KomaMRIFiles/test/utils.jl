function quantize_to_pulseq(seq::Sequence)
    blockDurationRaster = KomaMRIFiles.get_blockDurationRaster(seq)
    qseq = deepcopy(seq)
    for bi in 1:length(qseq)
        # RFs
        qseq.RF[1, bi].delay = quantize_time(qseq.RF[1, bi].delay, 1e-6) # quantize delay to us
        # Gradients
        for gi in 1:3
            qseq.GR[gi, bi].T     = quantize_time(qseq.GR[gi, bi].T, 1e-6) # quantize time to us
            qseq.GR[gi, bi].delay = quantize_time(qseq.GR[gi, bi].delay, 1e-6) # quantize delay to us
            qseq.GR[gi, bi].rise  = quantize_time(qseq.GR[gi, bi].rise, 1e-6) # quantize rise to us
            qseq.GR[gi, bi].fall  = quantize_time(qseq.GR[gi, bi].fall, 1e-6) # quantize fall to us
        end
        # ADCs
        qseq.ADC[bi].T     = quantize_time(qseq.ADC[bi].T, 1e-6) # quantize time to us
        qseq.ADC[bi].delay = quantize_time(qseq.ADC[bi].delay, 1e-6) # quantize delay to us
        # Durations
        max_event_duration = max(dur.(qseq.GR[:, bi])..., dur(qseq.RF[1, bi]), dur(qseq.ADC[bi]))
        block_duration = max(qseq.DUR[bi], max_event_duration)
        # Round to block duration raster (round up to ensure >= max_event_duration)
        qseq.DUR[bi] = ceil(block_duration / blockDurationRaster) * blockDurationRaster
    end
    return qseq
end

function round_trip_seq()
    seq = Sequence()
    seq += RF(1.0, 100e-6)
    seq += Grad(1.0, 100e-6, 10e-6)
    seq += ADC(100, 10e-6, 10e-6)
    return seq
end

quantize_time(t::Array, factor::Float64) = t
quantize_time(t::Real,  factor::Float64) = round(Int, t / factor) * factor