# ===========================================================================
# 6. Full Sequence Sampling
# ===========================================================================
#
# Assemble sampled blocks into the DiscreteSequence used by downstream sequence
# calculations.

# -- 6.1. Add external key times --------------------------------------------
sequence_boundary_sampling_times(seq) = (zero(dur(seq)), dur(seq))
motion_sampling_times(seq, ::NoMotion) = eltype(sequence_boundary_sampling_times(seq))[]
function motion_sampling_times(seq, motion)
    t = [zero(dur(seq)), dur(seq)]
    add_key_time_points!(t, motion)
    sort!(unique!(t))
    return t[2:(end - 1)]
end

block_global_event_times(T0, block, global_event_times) =
    [t - T0[block] for t in global_event_times if T0[block] <= t <= T0[block + 1]]

# -- 6.2. Assemble sampled blocks -------------------------------------------
function append_adc_start_padding!(out, values, first_t)
    first(values.ADC) || return out
    push!(out.t, first_t)
    foreach((dst, src) -> push!(dst, first(src)), table_columns(out)[1:6], table_columns(values)[1:6])
    push!(out.ADC, false)
    push!(out.Δt, 0.0)
    push!(out.excitation_bool, false)
    return out
end

function same_boundary_sample(out, values)
    return all(last(dst) == first(src) for (dst, src) in zip(table_columns(out), table_columns(values)))
end

function append_sampled_block!(out, values, t0)
    isempty(values.t) && return out
    first_t = t0 + first(values.t)
    first_row = firstindex(values.t)
    isempty(out.t) ? append_adc_start_padding!(out, values, first_t) : if first_t == last(out.t) && same_boundary_sample(out, values)
        first_row += 1
    else
        push!(out.Δt, first_t - last(out.t))
        push!(out.excitation_bool, false)
    end
    rows = first_row:lastindex(values.t)
    for t in view(values.t, rows)
        push!(out.t, t0 + t)
    end
    foreach((dst, src) -> append!(dst, view(src, rows)), table_columns(out), table_columns(values))
    append!(out.Δt, values.Δt)
    append!(out.excitation_bool, values.excitation_bool)
    return out
end

# -- 6.3. Sample the full sequence ------------------------------------------
function sample_sequence(seq; motion=NoMotion(), sampling_rule=MaxStepSizeRule(1e-3, 5e-5), freq_in_phase=false)
    out = DiscreteSequence()
    sizehint = max(8length(seq), 5sum(seq.ADC.N))
    foreach(x -> Base.sizehint!(x, sizehint), (out.t, table_columns(out)..., out.excitation_bool, out.Δt))
    T0 = get_block_start_times(seq)
    global_event_times = merge_sampling_times(sequence_boundary_sampling_times(seq), motion_sampling_times(seq, motion))
    for block in 1:length(seq)
        values = sample_sequence_block(seq, block; sampling_rule, motion_times=block_global_event_times(T0, block, global_event_times), freq_in_phase)
        append_sampled_block!(out, values, T0[block])
    end
    return out
end

# -- 6.4. Public discretization entrypoint ----------------------------------
"""
    seqd = discretize(seq::Sequence; sampling_rule=MaxStepSizeRule(1e-3, 5e-5))

This function returns a sampled Sequence struct with RF and gradient time refinements.

# Arguments
- `seq`: (`::Sequence`) sequence

# Keywords
- `sampling_rule`: controls how sequence sampling times are refined
- `freq_in_phase`: fold RF frequency modulation into the complex RF waveform

# Returns
- `seqd`: (`::DiscreteSequence`) DiscreteSequence struct
"""
function discretize(seq::Sequence; motion=NoMotion(), sampling_rule=MaxStepSizeRule(1e-3, 5e-5), freq_in_phase=false)
    return sample_sequence(seq; motion, sampling_rule, freq_in_phase)
end

# -- 6.5. Evaluate the sequence at prescribed times -------------------------
function evaluate_sequence_at(seq, t; freq_in_phase=false)
    flat_t = vec(t)
    out = DiscreteSequence(flat_t)
    isempty(flat_t) && return out
    T0 = get_block_start_times(seq)
    for block in 1:length(seq)
        idx = findall(t -> T0[block] <= t && (t < T0[block + 1] || block == length(seq)), flat_t)
        isempty(idx) && continue
        local_t = flat_t[idx] .- T0[block]
        p = sortperm(local_t; alg=MergeSort)
        rf = seq.RF[1, block]
        values = sample_sequence_block(
            block_event_waveforms(rf, seq.GR[1, block], seq.GR[2, block], seq.GR[3, block], seq.ADC[block]),
            local_t[p];
            rf_center=rf_center_time(rf),
            freq_in_phase,
        )
        foreach((dst, src) -> dst[idx[p]] .= src, table_columns(out), table_columns(values))
    end
    return out
end
