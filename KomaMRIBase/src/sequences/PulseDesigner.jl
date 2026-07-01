"""
    PulseDesigner

A module to define different pulse sequences.
"""
module PulseDesigner
using ..KomaMRIBase
using FFTW: fft, fftshift, ifftshift
using Interpolations: linear_interpolation

ramp_time_to_raster(t, raster) = max(ceil_to_raster(t, raster), raster)

slew_limited_ramp_samples(amplitude, max_slew, raster) =
    raster_samples(ceil_to_raster(abs(amplitude) / max_slew, raster), raster)

slew_limited_rise_time(amplitude; sys, max_slew) =
    ramp_time_to_raster(abs(amplitude) / max_slew, sys.GR_Δt)

normalize_flip_angle(signal, dwell, flip_angle) =
    flip_angle / (2π * γ * abs(sum(signal) * dwell)) .* signal

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
export ceil_to_raster, floor_to_raster, round_to_raster
export EPI, radial_base, EPI_example
end
