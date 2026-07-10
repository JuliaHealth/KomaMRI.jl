const SLR_COMPLEMENT_MARGIN = 1e-7
const SLR_SPECTRAL_FACTOR_MARGIN = 1e-6

"""
    seq = build_slr_pulse(flip_angle; kwargs...)

Return a `Sequence` with a Pulseq-style SLR pulse. See
`make_slr_pulse` for pulse design keywords.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Returns
- `seq`: Sequence containing the RF block and optional slice gradients.
"""
function build_slr_pulse(flip_angle; sys=Scanner(), kwargs...)
    rf, gz, gz_rephaser = make_slr_pulse(flip_angle; sys, kwargs...)
    seq = Sequence(sys)
    if gz === nothing
        addblock!(seq, rf)
        seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
        return seq
    end
    addblock!(seq, rf; z=gz)
    seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
    addblock!(seq; z=gz_rephaser)
    return seq
end

"""
    rf, gz, gzr = make_slr_pulse(flip_angle; duration=1e-3, sys=Scanner(), kwargs...)

Return a Shinnar-Le Roux RF event tuple matching Pulseq's `makeSLRpulse`
design paths. `gz` and `gzr` are `nothing` unless `slice_thickness` is
supplied.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Keywords
- `duration=1e-3`: RF pulse duration. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.
- `slice_thickness=nothing`: Slice thickness for slice-selective RF. [`m`]
- `time_bw_product=4.0`: Time-bandwidth product.
- `passband_ripple=0.01`: Passband ripple passed to the SLR design.
- `stopband_ripple=0.01`: Stopband ripple passed to the SLR design.
- `filter_type=:ms`: FIR design: `:ms`, `:pm`, `:min`, `:max`, or `:ls`.
- `freq_offset=0.0`: RF frequency offset. [`Hz`]
- `phase_offset=0.0`: RF phase offset. [`rad`]
- `delay=0.0`: RF delay before RF dead-time adjustment. [`s`]
- `dwell=sys.RF_Δt`: RF sample spacing. [`s`]
- `use=Excitation()`: RF use label controlling the SLR pulse type.
- `max_grad=nothing`: Slice-gradient amplitude limit override. Plain numbers use Pulseq units. [`Hz/m`]
- `max_slew=nothing`: Slice-gradient slew limit override. Plain numbers use Pulseq units. [`Hz/m/s`]

# Returns
- `rf`: RF event.
- `gz`: Slice-select gradient event.
- `gzr`: Slice rephaser gradient event.

# References
- Pauly, J., Le Roux, P., Nishimura, D. and Macovski, A. (1991). "Parameter
  Relations for the Shinnar-Le Roux Selective Excitation Pulse Design
  Algorithm". IEEE Transactions on Medical Imaging, 10(1), 53-65.
  https://doi.org/10.1109/42.75611
"""
function make_slr_pulse(flip_angle; duration=1e-3, sys=Scanner(),
    slice_thickness=nothing, time_bw_product=4.0, passband_ripple=0.01,
    stopband_ripple=0.01, filter_type=:ms, freq_offset=0.0, phase_offset=0.0,
    delay=0.0, dwell=sys.RF_Δt, use=Excitation(), max_grad=nothing,
    max_slew=nothing)
    flip_angle      = to_SI(flip_angle, SIUnitsDefault())
    duration        = to_SI(duration, SIUnitsDefault())
    slice_thickness = to_SI(slice_thickness, SIUnitsDefault())
    freq_offset     = to_SI(freq_offset, SIUnitsDefault())
    phase_offset    = to_SI(phase_offset, SIUnitsDefault())
    delay           = to_SI(delay, SIUnitsDefault())
    dwell           = to_SI(dwell, SIUnitsDefault())
    max_grad        = isnothing(max_grad) ? sys.Gmax : to_SI(max_grad, PulseqUnitsDefault())
    max_slew        = isnothing(max_slew) ? sys.Smax : to_SI(max_slew, PulseqUnitsDefault())
    duration > 0 || error("RF pulse duration must be positive.")
    dwell > 0 || error("RF dwell time must be positive.")
    time_bw_product > 0 || error("RF time-bandwidth product must be positive.")
    passband_ripple > 0 || error("Passband ripple must be positive.")
    stopband_ripple > 0 || error("Stopband ripple must be positive.")
    # MATLAB rounds positive half-sample durations away from zero.
    n = floor(Int, duration / dwell + 0.5)
    n > 0 || error("RF pulse duration is shorter than the RF raster.")
    waveform = slr_waveform(
        n, time_bw_product, use, flip_angle, Val(filter_type),
        passband_ripple, stopband_ripple,
    )
    waveform = normalize_flip_angle(waveform, dwell, flip_angle)
    center = rf_peak_center(waveform, dwell)
    rf_start_time = max(delay, sys.RF_dead_time)
    rf = RF(waveform, (n - 1) * dwell, freq_offset, rf_start_time + dwell / 2;
        center=center - dwell / 2, ϕ=phase_offset, use)
    slice_thickness === nothing && return rf, nothing, nothing
    slice_thickness > 0 || error("Slice thickness must be positive.")
    slice_area = time_bw_product / (γ * slice_thickness)
    gz, gz_rephaser = slice_select_gradient_events(duration, slice_area, rf; sys,
        rf_start_time, max_grad, max_slew)
    rf.delay = max(rf_start_time, gz.rise + gz.delay) + dwell / 2
    return rf, gz, gz_rephaser
end

function slr_waveform(n, time_bw_product, use, flip_angle, filter_type,
    passband_ripple, stopband_ripple)
    scale, passband_ripple, stopband_ripple, apply_inverse = slr_polynomial_spec(
        use, flip_angle, passband_ripple, stopband_ripple,
    )
    beta = slr_beta_polynomial(
        n, time_bw_product, filter_type, passband_ripple, stopband_ripple,
    )
    beta .*= scale
    return apply_inverse ? inverse_slr_transform(beta) : complex.(beta)
end

function slr_beta_polynomial(
    n, time_bw_product, ::Val{:ms}, passband_ripple, stopband_ripple,
)
    x = ((0:(n - 1)) .- n / 2) ./ (n / 2)
    return @. sinc(time_bw_product * x / 2) *
        (0.54 + 0.46cos(π * x)) * time_bw_product / n
end

function slr_beta_polynomial(
    n, time_bw_product, ::Val{:pm}, passband_ripple, stopband_ripple,
)
    width = slr_transition_measure(passband_ripple, stopband_ripple) /
        time_bw_product
    bands = [
        0.0,
        (1 - width) * time_bw_product / 2,
        (1 + width) * time_bw_product / 2,
        n / 2,
    ] ./ n
    return remez(
        n, bands, [1.0, 0.0]; weight=[1.0, passband_ripple / stopband_ripple],
    )
end

function slr_beta_polynomial(
    n, time_bw_product, ::Val{:min}, passband_ripple, stopband_ripple,
)
    return reverse(minimum_phase_polynomial(
        n, time_bw_product, passband_ripple, stopband_ripple,
    ))
end

function slr_beta_polynomial(
    n, time_bw_product, ::Val{:max}, passband_ripple, stopband_ripple,
)
    return minimum_phase_polynomial(
        n, time_bw_product, passband_ripple, stopband_ripple,
    )
end

function slr_beta_polynomial(
    n, time_bw_product, ::Val{:ls}, passband_ripple, stopband_ripple,
)
    iseven(n) || error("The :ls SLR filter requires an even number of samples.")
    width = slr_transition_measure(passband_ripple, stopband_ripple) /
        time_bw_product
    bands = [
        0.0,
        (1 - width) * time_bw_product / 2,
        (1 + width) * time_bw_product / 2,
        n / 2,
    ] ./ (n / 2)
    coefficients = least_squares_fir(n + 1, bands, [1.0, 0.0];
        weight=[1.0, passband_ripple / stopband_ripple])
    indices = [0:div(n, 2); -div(n, 2):-1]
    phase = cis.(2π .* indices ./ (2(n + 1)))
    shifted = ifft(fft(coefficients) .* phase)
    return real.(shifted[1:n])
end

slr_beta_polynomial(n, time_bw_product, ::Val{filter_type}, passband_ripple,
    stopband_ripple) where {filter_type} = error(
    "Unsupported SLR filter type `:$filter_type`; use :ms, :pm, :min, :max, or :ls.",
)

# Pauly et al., Table I: convert slice-profile ripple to B-polynomial ripple.
slr_polynomial_spec(::Excitation, flip_angle, passband_ripple, stopband_ripple) =
    flip_angle <= π / 6 ?
    (1.0, passband_ripple, stopband_ripple, false) :
    (sqrt(0.5), sqrt(passband_ripple / 2), stopband_ripple / sqrt(2), true)
slr_polynomial_spec(::Refocusing, flip_angle, passband_ripple, stopband_ripple) =
    (1.0, passband_ripple / 4, sqrt(stopband_ripple), true)
slr_polynomial_spec(::Inversion, flip_angle, passband_ripple, stopband_ripple) =
    (1.0, passband_ripple / 8, sqrt(stopband_ripple / 2), true)
slr_polynomial_spec(::Saturation, flip_angle, passband_ripple, stopband_ripple) =
    (sqrt(0.5), passband_ripple / 2, sqrt(stopband_ripple), true)
slr_polynomial_spec(::Preparation, flip_angle, passband_ripple, stopband_ripple) =
    (1.0, passband_ripple, stopband_ripple, false)
slr_polynomial_spec(::Other, flip_angle, passband_ripple, stopband_ripple) =
    (1.0, passband_ripple, stopband_ripple, false)
slr_polynomial_spec(::Undefined, flip_angle, passband_ripple, stopband_ripple) =
    (1.0, passband_ripple, stopband_ripple, false)

# Pauly et al., Eq. (21): empirical optimal-FIR transition measure.
function slr_transition_measure(passband_ripple, stopband_ripple)
    log_pass = log10(passband_ripple)
    log_stop = log10(stopband_ripple)
    return (5.309e-3log_pass^2 + 7.114e-2log_pass - 4.761e-1) * log_stop +
        (-2.66e-3log_pass^2 - 5.941e-1log_pass - 4.278e-1)
end

function minimum_phase_polynomial(n, time_bw_product, passband_ripple,
    stopband_ripple)
    adjusted_stopband = 0.5stopband_ripple^2
    width = 0.5slr_transition_measure(2passband_ripple, adjusted_stopband) /
        time_bw_product
    bands = [
        0.0,
        (1 - width) * time_bw_product / 2,
        (1 + width) * time_bw_product / 2,
        n / 2,
    ] ./ n
    linear_phase = remez(
        2n - 1, bands, [1.0, 0.0];
        weight=[1.0, 2passband_ripple / adjusted_stopband],
    )

    fft_length = 128nextpow(2, length(linear_phase))
    left_padding = cld(fft_length - length(linear_phase), 2)
    padded = zeros(Float64, fft_length)
    padded[(left_padding + 1):(left_padding + length(linear_phase))] = linear_phase
    power_response = fftshift(fft(ifftshift(complex.(padded))))
    positive_response = power_response .-
        (1 + SLR_SPECTRAL_FACTOR_MARGIN) * minimum(real, power_response)
    magnitude = sqrt.(abs.(positive_response))
    minimum_phase = minimum_phase_response(magnitude)
    time_response = ifft(ifftshift(conj.(minimum_phase)))
    return time_response[1:n]
end

function least_squares_fir(num_taps, bands, desired; weight)
    half_order = div(num_taps - 1, 2)
    band_pairs = ((bands[1], bands[2]), (bands[3], bands[4]))
    q = [sum(weight[j] * (
        last(band_pairs[j]) * sinc(k * last(band_pairs[j])) -
        first(band_pairs[j]) * sinc(k * first(band_pairs[j]))
    ) for j in eachindex(band_pairs)) for k in 0:(num_taps - 1)]
    gram = [q[abs(i - j) + 1] + q[i + j - 1]
        for i in 1:(half_order + 1), j in 1:(half_order + 1)]
    moment = [sum(weight[j] * desired[j] * (
        last(band_pairs[j]) * sinc(k * last(band_pairs[j])) -
        first(band_pairs[j]) * sinc(k * first(band_pairs[j]))
    ) for j in eachindex(band_pairs)) for k in 0:half_order]
    coefficients = cholesky(Symmetric(gram)) \ moment
    return [reverse(coefficients[2:end]); 2coefficients[1]; coefficients[2:end]]
end

# Pauly et al., Eqs. (15) and (16): remove one Cayley-Klein rotation at a time.
function inverse_slr_transform(beta)
    alpha = complementary_slr_polynomial(beta)
    beta = complex.(beta)
    waveform = zeros(ComplexF64, length(beta))
    for i in reverse(eachindex(beta))
        ratio = beta[i] / alpha[i]
        cosine = inv(sqrt(1 + abs2(ratio)))
        sine = conj(cosine * ratio)
        waveform[i] = 2atan(abs(sine), cosine) * cis(angle(sine))
        if i > 1
            previous_alpha = cosine .* alpha .+ sine .* beta
            previous_beta = -conj(sine) .* alpha .+ cosine .* beta
            alpha = previous_alpha[2:i]
            beta = previous_beta[1:(i - 1)]
        end
    end
    return waveform
end

# Pauly et al., Eqs. (14), (17), and (18): choose the minimum-energy
# complementary A polynomial satisfying |A|^2 + |B|^2 = 1.
function complementary_slr_polynomial(beta)
    n = length(beta)
    padded = zeros(ComplexF64, 16n)
    padded[1:n] = beta
    beta_response = fft(padded)
    peak = maximum(abs, beta_response)
    peak < 1 || (beta_response ./= peak + SLR_COMPLEMENT_MARGIN)
    alpha_magnitude = sqrt.(max.(0.0, 1 .- abs2.(beta_response)))
    alpha_response = minimum_phase_response(alpha_magnitude)
    alpha = fft(alpha_response) ./ length(padded)
    return reverse(alpha[1:n])
end

# Eq. (18), evaluated with the discrete analytic signal of log magnitude.
function minimum_phase_response(magnitude)
    n = length(magnitude)
    analytic_log_magnitude = fft(log.(abs.(magnitude)))
    analytic_log_magnitude[2:div(n, 2)] .*= 2
    analytic_log_magnitude[(div(n, 2) + 2):end] .= 0
    return exp.(ifft(analytic_log_magnitude))
end
