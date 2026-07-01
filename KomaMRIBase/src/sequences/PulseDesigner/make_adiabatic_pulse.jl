"""
    seq = build_adiabatic_pulse(type; kwargs...)

Return a `Sequence` with a Pulseq-style adiabatic RF pulse. See
`make_adiabatic_pulse` for pulse design keywords.

# Arguments
- `type`: Pulse type, `:hypsec` or `:wurst`.

# Returns
- `seq`: Sequence containing the RF block and optional slice gradients.
"""
function build_adiabatic_pulse(type; sys=Scanner(), kwargs...)
    rf, gz, gz_rephaser = make_adiabatic_pulse(type; sys, kwargs...)
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
    rf, gz, gzr = make_adiabatic_pulse(type; sys=Scanner(), kwargs...)

Return a Pulseq-style adiabatic RF event tuple. `gz` and `gzr` are `nothing`
unless `slice_thickness` is supplied.

# Arguments
- `type`: Pulse type, `:hypsec` or `:wurst`.

# Keywords
Common RF keywords:
- `duration=10e-3`: RF pulse duration. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.
- `freq_offset=0.0`: RF frequency offset. [`Hz`]
- `phase_offset=0.0`: RF phase offset. [`rad`]
- `adiabaticity=4.0`: Adiabaticity scaling.
- `delay=0.0`: RF delay before RF dead-time adjustment. [`s`]
- `dwell=sys.RF_Δt`: RF sample spacing. [`s`]
- `use=Undefined()`: RF use label.

Hypsec keywords:
- `beta=800.0`: AM waveform parameter.
- `mu=4.9`: FM sweep parameter.

WURST keywords:
- `n_fac=40`: AM exponent.
- `bandwidth=40000.0`: WURST bandwidth. [`Hz`]

Slice-selective keywords:
- `slice_thickness=nothing`: Slice thickness. [`m`]
- `max_grad=nothing`: Slice-gradient amplitude limit override. Plain numbers use Pulseq units. [`Hz/m`]
- `max_slew=nothing`: Slice-gradient slew limit override. Plain numbers use Pulseq units. [`Hz/m/s`]

# Returns
- `rf`: RF event.
- `gz`: Slice-select gradient event.
- `gzr`: Slice rephaser gradient event.

# References
- `:hypsec` from Baum, J., Tycko, R. and Pines, A. (1985). "Broadband and
  adiabatic inversion of a two-level system by phase-modulated pulses".
  Phys. Rev. A, 32:3435-3447.
- `:wurst` from Kupce, E. and Freeman, R. (1995). "Stretched Adiabatic Pulses
  for Broadband Spin Inversion". J. Magn. Reson. Ser. A, 117:246-256.
"""
make_adiabatic_pulse(type; kwargs...) =
    error("Adiabatic pulse type must be a Symbol, for example :hypsec or :wurst.")
make_adiabatic_pulse(type::Symbol; kwargs...) = make_adiabatic_pulse(Val(type); kwargs...)
make_adiabatic_pulse(type::Val; kwargs...) = _make_adiabatic_pulse(type; kwargs...)
_make_adiabatic_pulse(::Val{type}; kwargs...) where {type} =
    error("Unsupported adiabatic pulse type `:$type`; use :hypsec or :wurst.")

function _make_adiabatic_pulse(type::Union{Val{:hypsec},Val{:wurst}};
    duration=10e-3, sys=Scanner(), slice_thickness=nothing,
    freq_offset=0.0, phase_offset=0.0, beta=800.0, mu=4.9, n_fac=40,
    bandwidth=40000.0, adiabaticity=4.0, delay=0.0, dwell=sys.RF_Δt,
    use=Undefined(), max_grad=nothing, max_slew=nothing)
    duration        = to_SI(duration, SIUnitsDefault())
    slice_thickness = to_SI(slice_thickness, SIUnitsDefault())
    freq_offset     = to_SI(freq_offset, SIUnitsDefault())
    phase_offset    = to_SI(phase_offset, SIUnitsDefault())
    bandwidth       = to_SI(bandwidth, SIUnitsDefault())
    delay           = to_SI(delay, SIUnitsDefault())
    dwell           = to_SI(dwell, SIUnitsDefault())
    max_grad        = isnothing(max_grad) ? sys.Gmax : to_SI(max_grad, PulseqUnitsDefault())
    max_slew        = isnothing(max_slew) ? sys.Smax : to_SI(max_slew, PulseqUnitsDefault())
    duration > 0 || error("RF pulse duration must be positive.")
    dwell > 0 || error("RF dwell time must be positive.")
    nraw = round(Int, duration / dwell + eps())
    n = 4 * div(nraw, 4)
    n > 0 || error("RF pulse duration is shorter than four RF raster samples.")
    B₁, Δω = adiabatic_modulation(type, n, duration; beta, mu, n_fac, bandwidth)
    signal = adiabatic_signal(B₁, Δω, dwell, adiabaticity)
    if n != nraw
        npad = nraw - n
        signal = [zeros(ComplexF64, npad - div(npad, 2)); signal; zeros(ComplexF64, div(npad, 2))]
    end
    peak_center = rf_peak_center(signal, dwell)
    rf_start_time = max(delay, sys.RF_dead_time)
    rf = RF(signal ./ γ, (nraw - 1) * dwell, freq_offset, rf_start_time + dwell / 2;
        center=peak_center - dwell / 2, ϕ=phase_offset, use)
    slice_thickness === nothing && return rf, nothing, nothing
    slice_thickness > 0 || error("Slice thickness must be positive.")
    slice_bandwidth = adiabatic_slice_bandwidth(
        type, signal, peak_center, dwell, duration, bandwidth,
        freq_offset, phase_offset, Δω,
    )
    slice_area = slice_bandwidth * duration / (γ * slice_thickness)
    gz, gz_rephaser = slice_select_gradient_events(duration, slice_area, rf; sys,
        rf_start_time, max_grad, max_slew)
    rf.delay = max(rf_start_time, gz.rise + gz.delay) + dwell / 2
    return rf, gz, gz_rephaser
end

function adiabatic_modulation(::Val{:hypsec}, N, T; beta, mu, kwargs...)
    # Baum/Tycko/Pines HS pulse: B₁(t) = sech(βt), Δω(t) = -μβ tanh(βt).
    β, μ = beta, mu
    t = ((-div(N, 2)):(div(N, 2) - 1)) ./ N .* T
    B₁ = @. 1 / cosh(β * t)
    Δω = @. -μ * β * tanh(β * t)
    return B₁, Δω
end

function adiabatic_modulation(::Val{:wurst}, N, T; n_fac, bandwidth, kwargs...)
    # Kupce/Freeman WURST-n pulse: B₁(t) = 1 - |cos(πt/T)|ⁿ, linear Δω(t).
    n = n_fac
    t = (0:(N - 1)) .* T ./ N
    B₁ = @. 1 - abs(cospi(t / T))^n
    Δω = collect(range(-bandwidth / 2, bandwidth / 2; length=N)) .* 2π
    return B₁, Δω
end

function adiabatic_signal(B₁, Δω, Δt, adiabaticity)
    # Scale B₁ from the adiabatic condition at the sweep center:
    # (2π * B₁_Hz)^2 / abs(dΔω/dt) = adiabaticity.
    ϕ = [zero(eltype(Δω)); cumtrapz(fill(Δt, length(Δω) - 1), Δω)]
    i = argmin(abs.(Δω))
    if iszero(Δω[i])
        ϕ0 = ϕ[i]
        B₁_0 = B₁[i]
        dΔω_dt = abs(Δω[i + 1] - Δω[i - 1]) / (2Δt)
    else
        step = Δω[i] * Δω[i + 1] < 0 ? 1 : -1
        ϕ0 = (ϕ[i] * Δω[i + step] - ϕ[i + step] * Δω[i]) / (Δω[i + step] - Δω[i])
        B₁_0 = (B₁[i] * Δω[i + step] - B₁[i + step] * Δω[i]) / (Δω[i + step] - Δω[i])
        dΔω_dt = abs(Δω[i] - Δω[i + step]) / Δt
    end
    B₁_Hz = sqrt(dΔω_dt * adiabaticity) / (2π * B₁_0)
    return B₁_Hz .* B₁ .* cis.(ϕ .- ϕ0)
end

function rf_peak_center(signal, dwell)
    # MATLAB Pulseq calcRfCenter uses the peak plateau center; Koma's default
    # RF center is amplitude-weighted, so choose the Pulseq center explicitly.
    peak = maximum(abs, signal)
    idx = findall(abs.(signal) .>= peak * 0.99999)
    t = ((1:length(signal)) .- 0.5) .* dwell
    return (t[first(idx)] + t[last(idx)]) / 2
end

adiabatic_slice_bandwidth(::Val{:wurst}, signal, center, dwell, duration, bandwidth,
    freq_offset, phase_offset, Δω) = bandwidth

function adiabatic_slice_bandwidth(::Val{:hypsec}, signal, center, dwell, duration,
    bandwidth, freq_offset, phase_offset, Δω)
    # MATLAB Pulseq uses calcRfBandwidth(rf, 0.1) for hypsec slice gradients.
    return rf_bandwidth(signal, center, dwell, freq_offset, phase_offset, 0.1)
end

function rf_bandwidth(signal, center, dwell, freq_offset, phase_offset, cutoff)
    df = 10.0
    dt = 1e-6
    rf_t = ((1:length(signal)) .- 0.5) .* dwell
    values = signal .* cis.(phase_offset .+ 2π .* freq_offset .* rf_t)
    centered_t = rf_t .- center
    n = round(1 / df / dt)
    sample_t = (-floor(n / 2):(ceil(n / 2) - 1)) .* dt
    sampled = linear_interpolation(centered_t, values; extrapolation_bc=0.0).(sample_t)
    spectrum = abs.(fftshift(fft(ifftshift(sampled))))
    f = (-floor(n / 2):(ceil(n / 2) - 1)) .* df
    return rf_flank(reverse(f), reverse(spectrum), cutoff) - rf_flank(f, spectrum, cutoff)
end

function rf_flank(f, spectrum, cutoff)
    s = spectrum ./ maximum(spectrum) .- cutoff
    i = findfirst(>(0), s)
    i === nothing && error("Could not determine RF bandwidth.")
    i == firstindex(s) && return f[i]
    s0, s1 = s[i - 1], s[i]
    return (s1 * f[i - 1] - s0 * f[i]) / (s1 - s0)
end
