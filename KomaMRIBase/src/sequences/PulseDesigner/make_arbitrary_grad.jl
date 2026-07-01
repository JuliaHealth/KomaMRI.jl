"""
    seq = build_arbitrary_grad(axis::Symbol, waveform; kwargs...)

Return a one-block `Sequence` with a Pulseq-style arbitrary gradient on `axis`.
See `make_arbitrary_grad` for gradient design keywords.

# Arguments
- `axis::Symbol`: Gradient axis, one of `:x`, `:y`, or `:z`.
- `waveform`: Gradient waveform samples. Plain numbers use Pulseq units. [`Hz/m`]

# Returns
- `seq`: Sequence containing the gradient block.
"""
function build_arbitrary_grad(axis::Symbol, waveform; sys=Scanner(), kwargs...)
    seq = Sequence(sys)
    axis in (:x, :y, :z) || error("Gradient axis must be :x, :y, or :z.")
    grad = make_arbitrary_grad(waveform; sys, kwargs...)
    addblock!(seq; (; axis => grad)...)
    return seq
end

"""
    grad = make_arbitrary_grad(waveform; kwargs...)

Return a Pulseq-style arbitrary gradient event.

# Arguments
- `waveform`: Gradient waveform samples. Plain numbers use Pulseq units. [`Hz/m`]

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.
- `oversampling=false`: Use Pulseq oversampled arbitrary-gradient timing.
- `max_grad=nothing`: Gradient amplitude limit override. Plain numbers use Pulseq units. [`Hz/m`]
- `max_slew=nothing`: Slew limit override. Plain numbers use Pulseq units. [`Hz/m/s`]
- `delay=0.0`: Gradient delay. [`s`]
- `first=nothing`: Gradient amplitude before the first sample. Plain numbers use Pulseq units. [`Hz/m`]
- `last=nothing`: Gradient amplitude after the last sample. Plain numbers use Pulseq units. [`Hz/m`]

# Returns
- `grad`: Arbitrary gradient event.
"""
function make_arbitrary_grad(waveform; sys=Scanner(), oversampling=false,
    max_grad=nothing, max_slew=nothing, delay=0.0, first=nothing, last=nothing)
    waveform = to_SI.(collect(waveform), Ref(PulseqUnitsDefault()))
    max_grad = isnothing(max_grad) ? sys.Gmax : to_SI(max_grad, PulseqUnitsDefault())
    max_slew = isnothing(max_slew) ? sys.Smax : to_SI(max_slew, PulseqUnitsDefault())
    delay    = to_SI(delay, SIUnitsDefault())
    first    = to_SI(first, PulseqUnitsDefault())
    last     = to_SI(last, PulseqUnitsDefault())
    isempty(waveform) && error("Gradient waveform must not be empty.")
    if oversampling && iseven(length(waveform))
        error("Oversampled gradient waveforms must have odd length.")
    end
    if first === nothing || last === nothing
        length(waveform) > 1 || error("Supply first and last for one-sample gradients.")
    end
    if first === nothing
        first = oversampling ? 2 * waveform[1] - waveform[2] :
            (3 * waveform[1] - waveform[2]) / 2
    end
    if last === nothing
        last = oversampling ? 2 * waveform[end] - waveform[end - 1] :
            (3 * waveform[end] - waveform[end - 1]) / 2
    end
    duration = (length(waveform) - 1) * (oversampling ? sys.GR_Δt / 2 : sys.GR_Δt)
    grad = Grad(waveform, duration, sys.GR_Δt / 2, sys.GR_Δt / 2, delay, first, last)
    check_hw_limits(grad; max_grad, max_slew)
    return grad
end
