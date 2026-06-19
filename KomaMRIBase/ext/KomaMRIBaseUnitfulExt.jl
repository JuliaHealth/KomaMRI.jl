module KomaMRIBaseUnitfulExt

using Unitful
import Unitful: Frequency, Length, Time
import KomaMRIBase: PulseDesigner, Scanner, Undefined, γ
import KomaMRIBase:
    scanner_field_strength,
    scanner_gradient_amplitude,
    scanner_rf_amplitude,
    scanner_slew_rate,
    scanner_time

# Canonical SI and Pulseq-compatible units accepted by the wrappers.
const s                 = u"s"
const m                 = u"m"
const rad               = u"rad"
const deg               = u"deg"
const Hz                = u"Hz"
const tesla             = u"T"
const T_per_m           = u"T/m"
const Ts_per_m          = u"T*s/m"
const T_per_m_per_s     = u"T/m/s"
const Hz_per_m          = u"Hz/m"
const Hzs_per_m         = u"Hz*s/m"
const Hz_per_m_per_s    = u"Hz/m/s"

# Unitful quantity categories used for dispatch instead of runtime dimension checks.
const Radian                     = Quantity{T,NoDims,typeof(rad)} where {T}
const Degree                     = Quantity{T,NoDims,typeof(deg)} where {T}
const UnitfulAngle               = Union{Radian,Degree}
const MagneticFluxDensity        = Quantity{T,dimension(tesla),U} where {T,U}
const GradientAmplitude          = Quantity{T,dimension(T_per_m),U} where {T,U}
const PulseqGradientAmplitude    = Quantity{T,dimension(Hz_per_m),U} where {T,U}
const GradientArea               = Quantity{T,dimension(Ts_per_m),U} where {T,U}
const PulseqGradientArea         = Quantity{T,dimension(Hzs_per_m),U} where {T,U}
const SlewRate                   = Quantity{T,dimension(T_per_m_per_s),U} where {T,U}
const PulseqSlewRate             = Quantity{T,dimension(Hz_per_m_per_s),U} where {T,U}
const AnyGradientAmplitude       = Union{GradientAmplitude,PulseqGradientAmplitude}
const AnyGradientArea            = Union{GradientArea,PulseqGradientArea}
const RFAmplitude                = Union{MagneticFluxDensity,Frequency}
const UnitfulTimes               = AbstractVector{<:Time}
const UnitfulGradientAmplitudes  = AbstractVector{<:AnyGradientAmplitude}

# Strip Unitful values to Koma SI numbers; plain reals already mean SI.
seconds(::Nothing) = nothing
seconds(x::Real) = x
seconds(x::Time) = ustrip(Float64, s, x)

meters(::Nothing) = nothing
meters(x::Real) = x
meters(x::Length) = ustrip(Float64, m, x)

hertz(::Nothing) = nothing
hertz(x::Real) = x
hertz(x::Frequency) = ustrip(Float64, Hz, x)

radians(x::Real) = x
radians(x::UnitfulAngle) = ustrip(Float64, rad, x)

field_strength(x::Real) = x
field_strength(x::MagneticFluxDensity) = ustrip(Float64, tesla, x)

rf_amplitude(x::Real) = x
rf_amplitude(x::MagneticFluxDensity) = ustrip(Float64, tesla, x)
rf_amplitude(x::Frequency) = ustrip(Float64, Hz, x) / γ

gradient(::Nothing) = nothing
gradient(x::Real) = x
gradient(x::GradientAmplitude) = ustrip(Float64, T_per_m, x)
gradient(x::PulseqGradientAmplitude) = ustrip(Float64, Hz_per_m, x) / γ

gradient_area(x::GradientArea) = ustrip(Float64, Ts_per_m, x)
gradient_area(x::PulseqGradientArea) = ustrip(Float64, Hzs_per_m, x) / γ

slew_rate(x::Real) = x
slew_rate(x::SlewRate) = ustrip(Float64, T_per_m_per_s, x)
slew_rate(x::PulseqSlewRate) = ustrip(Float64, Hz_per_m_per_s, x) / γ

# Let Scanner's constructor convert unitful keyword values without owning Unitful.
scanner_field_strength(x::MagneticFluxDensity) = field_strength(x)
scanner_rf_amplitude(x::RFAmplitude) = rf_amplitude(x)
scanner_gradient_amplitude(x::AnyGradientAmplitude) = gradient(x)
scanner_slew_rate(x::Union{SlewRate,PulseqSlewRate}) = slew_rate(x)
scanner_time(x::Time) = seconds(x)

# Convert PulseDesigner inputs from Unitful to SI before calling core make_* methods.
# The build_* methods are unit-compatible because they call make_* internally.

PulseDesigner.si_time(x::Time) = seconds(x)
PulseDesigner.si_gradient(x::AnyGradientAmplitude) = gradient(x)
PulseDesigner.si_gradient_area(x::AnyGradientArea) = gradient_area(x)
PulseDesigner.si_slew_rate(x::Union{SlewRate,PulseqSlewRate}) = slew_rate(x)

function PulseDesigner.make_extended_trapezoid(
    times::UnitfulTimes, amplitudes::UnitfulGradientAmplitudes;
    sys=Scanner(), max_grad=sys.Gmax, max_slew=sys.Smax, kwargs...,
)
    times      = seconds.(times)
    amplitudes = gradient.(amplitudes)
    max_grad   = gradient(max_grad)
    max_slew   = slew_rate(max_slew)
    return PulseDesigner.make_extended_trapezoid(
        times, amplitudes; sys, max_grad, max_slew, kwargs...,
    )
end

function PulseDesigner.make_extended_trapezoid_area(
    grad_start::AnyGradientAmplitude, grad_end::AnyGradientAmplitude, area::AnyGradientArea;
    sys=Scanner(),
)
    grad_start = gradient(grad_start)
    grad_end   = gradient(grad_end)
    area       = gradient_area(area)
    return PulseDesigner.make_extended_trapezoid_area(grad_start, grad_end, area; sys)
end

function PulseDesigner.make_arbitrary_grad(waveform::UnitfulGradientAmplitudes;
    sys=Scanner(), max_grad=sys.Gmax, max_slew=sys.Smax, delay=0.0, first=nothing,
    last=nothing, kwargs...)
    waveform = gradient.(waveform)
    max_grad = gradient(max_grad)
    max_slew = slew_rate(max_slew)
    delay    = seconds(delay)
    first    = gradient(first)
    last     = gradient(last)
    return PulseDesigner.make_arbitrary_grad(
        waveform; sys, max_grad, max_slew, delay, first, last, kwargs...,
    )
end

function PulseDesigner.make_block_pulse(flip_angle::UnitfulAngle; duration=nothing,
    bandwidth=nothing, time_bw_product=0.0, sys=Scanner(), freq_offset=0.0,
    phase_offset=0.0, delay=0.0, use=Undefined())
    flip_angle   = radians(flip_angle)
    duration     = seconds(duration)
    bandwidth    = hertz(bandwidth)
    freq_offset  = hertz(freq_offset)
    phase_offset = radians(phase_offset)
    delay        = seconds(delay)
    return PulseDesigner.make_block_pulse(
        flip_angle; duration, bandwidth, time_bw_product, sys, freq_offset, phase_offset,
        delay, use,
    )
end

function PulseDesigner.make_sinc_pulse(flip_angle::UnitfulAngle; duration, sys=Scanner(),
    slice_thickness=nothing, freq_offset=0.0, phase_offset=0.0, time_bw_product=4.0,
    apodization=0.0, center_pos=0.5, delay=0.0, dwell=sys.RF_Δt, use=Undefined(),
    max_grad=sys.Gmax, max_slew=sys.Smax)
    flip_angle      = radians(flip_angle)
    duration        = seconds(duration)
    slice_thickness = meters(slice_thickness)
    freq_offset     = hertz(freq_offset)
    phase_offset    = radians(phase_offset)
    delay           = seconds(delay)
    dwell           = seconds(dwell)
    max_grad        = gradient(max_grad)
    max_slew        = slew_rate(max_slew)
    return PulseDesigner.make_sinc_pulse(
        flip_angle; duration, sys, slice_thickness, freq_offset, phase_offset,
        time_bw_product, apodization, center_pos, delay, dwell, use, max_grad, max_slew,
    )
end

function PulseDesigner.make_arbitrary_rf(signal, flip_angle::UnitfulAngle; sys=Scanner(),
    slice_thickness=nothing, bandwidth=nothing, time_bw_product=0.0, freq_offset=0.0,
    phase_offset=0.0, delay=0.0, dwell=sys.RF_Δt, center=nothing, use=Undefined(),
    max_grad=sys.Gmax, max_slew=sys.Smax)
    flip_angle      = radians(flip_angle)
    slice_thickness = meters(slice_thickness)
    bandwidth       = hertz(bandwidth)
    freq_offset     = hertz(freq_offset)
    phase_offset    = radians(phase_offset)
    delay           = seconds(delay)
    dwell           = seconds(dwell)
    center          = seconds(center)
    max_grad        = gradient(max_grad)
    max_slew        = slew_rate(max_slew)
    return PulseDesigner.make_arbitrary_rf(
        signal, flip_angle; sys, slice_thickness, bandwidth, time_bw_product,
        freq_offset, phase_offset, delay, dwell, center, use, max_grad, max_slew,
    )
end

function PulseDesigner.make_gauss_pulse(flip_angle::UnitfulAngle; duration, sys=Scanner(),
    slice_thickness=nothing, bandwidth=nothing, time_bw_product=3.0, freq_offset=0.0,
    phase_offset=0.0, apodization=0.0, center_pos=0.5, delay=0.0, dwell=sys.RF_Δt,
    use=Undefined(), max_grad=sys.Gmax, max_slew=sys.Smax)
    flip_angle      = radians(flip_angle)
    duration        = seconds(duration)
    slice_thickness = meters(slice_thickness)
    bandwidth       = hertz(bandwidth)
    freq_offset     = hertz(freq_offset)
    phase_offset    = radians(phase_offset)
    delay           = seconds(delay)
    dwell           = seconds(dwell)
    max_grad        = gradient(max_grad)
    max_slew        = slew_rate(max_slew)
    return PulseDesigner.make_gauss_pulse(
        flip_angle; duration, sys, slice_thickness, bandwidth, time_bw_product,
        freq_offset, phase_offset, apodization, center_pos, delay, dwell, use,
        max_grad, max_slew,
    )
end

function PulseDesigner.make_adiabatic_pulse(
    type::Union{Val{:hypsec},Val{:wurst}}; duration=10e-3, sys=Scanner(),
    slice_thickness=nothing, freq_offset=0.0, phase_offset=0.0, beta=800.0,
    mu=4.9, n_fac=40, bandwidth=40000.0, adiabaticity=4.0, delay=0.0,
    dwell=sys.RF_Δt, use=Undefined(), max_grad=sys.Gmax, max_slew=sys.Smax,
)
    duration        = seconds(duration)
    slice_thickness = meters(slice_thickness)
    freq_offset     = hertz(freq_offset)
    phase_offset    = radians(phase_offset)
    bandwidth       = hertz(bandwidth)
    delay           = seconds(delay)
    dwell           = seconds(dwell)
    max_grad        = gradient(max_grad)
    max_slew        = slew_rate(max_slew)
    return PulseDesigner._make_adiabatic_pulse(
        type; duration, sys, slice_thickness, freq_offset, phase_offset, beta, mu,
        n_fac, bandwidth, adiabaticity, delay, dwell, use, max_grad, max_slew,
    )
end

function PulseDesigner.make_rotation(phi::UnitfulAngle)
    return PulseDesigner.make_rotation(radians(phi))
end

function PulseDesigner.make_rotation(phi::UnitfulAngle, theta::UnitfulAngle)
    return PulseDesigner.make_rotation(radians(phi), radians(theta))
end

function PulseDesigner.make_rotation(phi::UnitfulAngle, theta::Real)
    return PulseDesigner.make_rotation(radians(phi), theta)
end

function PulseDesigner.make_rotation(phi::Real, theta::UnitfulAngle)
    return PulseDesigner.make_rotation(phi, radians(theta))
end

function PulseDesigner.make_rotation(axis::AbstractVector, angle::UnitfulAngle)
    return PulseDesigner.make_rotation(axis, radians(angle))
end

function PulseDesigner.make_trigger(channel::Symbol; delay=0.0, duration=0.0, sys=Scanner())
    delay    = seconds(delay)
    duration = seconds(duration)
    return PulseDesigner.trigger_event(PulseDesigner.trigger_channel(channel), delay, duration, sys)
end

function PulseDesigner.make_digital_output_pulse(
    channel::Symbol; delay=0.0, duration=0.0, sys=Scanner(),
)
    delay    = seconds(delay)
    duration = seconds(duration)
    return PulseDesigner.digital_output_event(
        PulseDesigner.digital_output_channel(channel), delay, duration, sys,
    )
end

function PulseDesigner.make_delay(delay::Time)
    delay = seconds(delay)
    return PulseDesigner.make_delay(delay)
end

function PulseDesigner.make_adc(num_samples, dwell::Time; duration=nothing, delay=0.0,
    sys=Scanner(), freq_offset=0.0, phase_offset=0.0)
    dwell        = seconds(dwell)
    duration     = seconds(duration)
    delay        = seconds(delay)
    freq_offset  = hertz(freq_offset)
    phase_offset = radians(phase_offset)
    return PulseDesigner.make_adc(
        num_samples, dwell; duration, delay, sys, freq_offset, phase_offset,
    )
end

# Keyword arguments do not dispatch, so use Integer to avoid replacing the core fallback.
function PulseDesigner.make_adc(num_samples::Integer; dwell=nothing, duration=nothing, delay=0.0,
    sys=Scanner(), freq_offset=0.0, phase_offset=0.0)
    dwell        = seconds(dwell)
    duration     = seconds(duration)
    delay        = seconds(delay)
    freq_offset  = hertz(freq_offset)
    phase_offset = radians(phase_offset)
    return PulseDesigner.make_adc(
        num_samples, dwell; duration, delay, sys, freq_offset, phase_offset,
    )
end

end
