"""
    PulseDesigner

A module to define different pulse sequences.
"""
module PulseDesigner
using ..KomaMRIBase
using FFTW: fft, fftshift, ifftshift
using Interpolations: linear_interpolation

function ceil_to_raster(t, raster)
    rounded = round(t / raster) * raster
    isapprox(t, rounded; rtol=0, atol=KomaMRIBase.PULSEQ_TIME_TOL) && return rounded
    return ceil(t / raster - KomaMRIBase.PULSEQ_TIME_TOL / raster) * raster
end

round_to_raster(t, raster) = round(Int, t / raster) * raster

raster_samples(t, raster) = round(Int, t / raster)

ramp_time_to_raster(t, raster) = max(ceil_to_raster(t, raster), raster)

slew_limited_ramp_samples(amplitude, max_slew, raster) =
    raster_samples(ceil_to_raster(abs(amplitude) / max_slew, raster), raster)

slew_limited_rise_time(amplitude; sys, max_slew) =
    ramp_time_to_raster(abs(amplitude) / max_slew, sys.GR_Δt)

si_time(::Nothing) = nothing
si_time(x) = x
si_gradient(::Nothing) = nothing
si_gradient(x) = x
si_gradient_area(::Nothing) = nothing
si_gradient_area(x) = x
si_slew_rate(x) = x

normalize_flip_angle(signal, dwell, flip_angle) =
    flip_angle / (2π * γ * abs(sum(signal) * dwell)) .* signal

function block_duration_to_fit_rf_ringdown(rf, sys, events...)
    rf_end_with_ringdown = KomaMRIBase._pulseq_duration(rf, sys.RF_Δt) +
        sys.RF_ring_down_time
    block_end = max(rf_end_with_ringdown, maximum(dur, events; init=0.0))
    return ceil_to_raster(block_end, sys.DUR_Δt)
end

include("PulseDesigner/make_trapezoid.jl")
include("PulseDesigner/make_arbitrary_grad.jl")
include("PulseDesigner/make_extended_trapezoid.jl")
include("PulseDesigner/make_extended_trapezoid_area.jl")
include("PulseDesigner/make_block_pulse.jl")
include("PulseDesigner/make_sinc_pulse.jl")
include("PulseDesigner/make_arbitrary_rf.jl")
include("PulseDesigner/make_gauss_pulse.jl")
include("PulseDesigner/make_adiabatic_pulse.jl")
include("PulseDesigner/make_label.jl")
include("PulseDesigner/make_rotation.jl")
include("PulseDesigner/make_trigger.jl")
include("PulseDesigner/make_digital_output_pulse.jl")
include("PulseDesigner/make_delay.jl")
include("PulseDesigner/make_adc.jl")
include("PulseDesigner/deprecated.jl")

export make_trapezoid, build_trapezoid
export make_arbitrary_grad, build_arbitrary_grad
export make_extended_trapezoid, build_extended_trapezoid
export make_extended_trapezoid_area, build_extended_trapezoid_area
export make_block_pulse, build_block_pulse
export make_sinc_pulse, build_sinc_pulse
export make_arbitrary_rf, build_arbitrary_rf
export make_gauss_pulse, build_gauss_pulse
export make_adiabatic_pulse, build_adiabatic_pulse
export make_label, build_label
export make_rotation, build_rotation
export make_trigger, build_trigger
export make_digital_output_pulse, build_digital_output_pulse
export make_delay, build_delay
export make_adc, build_adc
export EPI, radial_base, EPI_example
end
