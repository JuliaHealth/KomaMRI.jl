module KomaMRIBaseUnitfulExt

using Unitful
import Unitful: Frequency, Time
import KomaMRIBase: γ, to_SI
import KomaMRIBase:
    ceil_to_raster,
    floor_to_raster,
    round_to_raster,
    raster_samples,
    scanner_field_strength,
    scanner_gradient_amplitude,
    scanner_rf_amplitude,
    scanner_slew_rate,
    scanner_time

const s              = u"s"
const rad            = u"rad"
const deg            = u"deg"
const tesla          = u"T"
const T_per_m        = u"T/m"
const T_per_m_per_s  = u"T/m/s"
const Hz_per_m       = u"Hz/m"
const Hzs_per_m      = u"Hz*s/m"
const Hz_per_m_per_s = u"Hz/m/s"

const Radian                  = Quantity{T,NoDims,typeof(rad)} where {T}
const Degree                  = Quantity{T,NoDims,typeof(deg)} where {T}
const UnitfulAngle            = Union{Radian,Degree}
const MagneticFluxDensity     = Quantity{T,dimension(tesla),U} where {T,U}
const GradientAmplitude       = Quantity{T,dimension(T_per_m),U} where {T,U}
const PulseqGradientAmplitude = Quantity{T,dimension(Hz_per_m),U} where {T,U}
const PulseqGradientArea      = Quantity{T,dimension(Hzs_per_m),U} where {T,U}
const SlewRate                = Quantity{T,dimension(T_per_m_per_s),U} where {T,U}
const PulseqSlewRate          = Quantity{T,dimension(Hz_per_m_per_s),U} where {T,U}
const PulseqQuantities        = Union{
    PulseqGradientAmplitude,
    PulseqGradientArea,
    PulseqSlewRate,
}

to_SI(x::Quantity) = float(ustrip(upreferred(x)))
to_SI(x::UnitfulAngle) = float(ustrip(rad, x))
to_SI(x::Quantity{T,NoDims,U}) where {T,U} =
    throw(DimensionMismatch("Expected angle in rad or deg, got $(unit(x))."))
to_SI(x::PulseqQuantities) = float(ustrip(upreferred(x))) / γ

scanner_field_strength(x::MagneticFluxDensity) = to_SI(x)
scanner_rf_amplitude(x::MagneticFluxDensity) = to_SI(x)
scanner_rf_amplitude(x::Frequency) = to_SI(x) / γ
scanner_gradient_amplitude(x::GradientAmplitude) = to_SI(x)
scanner_slew_rate(x::SlewRate) = to_SI(x)
scanner_time(x::Time) = to_SI(x)

ceil_to_raster(t::Time, raster::Time) =
    ceil_to_raster(to_SI(t), to_SI(raster)) * s
floor_to_raster(t::Time, raster::Time) =
    floor_to_raster(to_SI(t), to_SI(raster)) * s
round_to_raster(t::Time, raster::Time) =
    round_to_raster(to_SI(t), to_SI(raster)) * s
raster_samples(t::Time, raster::Time) =
    raster_samples(to_SI(t), to_SI(raster))

end
