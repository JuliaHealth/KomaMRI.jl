module PulseDesignerUnitfulExt

using KomaMRIBase, Unitful

Angle{T} = Union{Quantity{T, NoDims, typeof(u"rad")}, Quantity{T, NoDims, typeof(u"°")}} where T

function PulseDesigner.block_pulse(flip_angle::Angle, duration::Unitful.Time;
    phase_offset::Angle=0u"°", freq_offset::Unitful.Frequency=0u"Hz", delay::Unitful.Time=0u"s", sys=Scanner())
    flip_angle, duration, phase_offset, freq_offset, delay = upreferred.((flip_angle, duration, phase_offset, freq_offset, delay)) .|> ustrip
    return PulseDesigner.block_pulse(FlipAngle(flip_angle), Duration(duration); phase_offset, freq_offset, delay, sys)
end

function PulseDesigner.block_pulse(flip_angle::Angle, bandwidth::Unitful.Frequency;
    phase_offset::Angle=0u"°", freq_offset::Unitful.Frequency=0u"Hz", delay::Unitful.Time=0u"s", sys=Scanner())
    flip_angle, bandwidth, phase_offset, freq_offset, delay = upreferred.((flip_angle, bandwidth, phase_offset, freq_offset, delay)) .|> ustrip
    return PulseDesigner.block_pulse(FlipAngle(flip_angle), Bandwidth(bandwidth); phase_offset, freq_offset, delay, sys)
end

function PulseDesigner.block_pulse(flip_angle::Angle, bandwidth::Unitful.Frequency, time_bw_product::DimensionlessQuantity;
    phase_offset::Angle=0u"°", freq_offset::Unitful.Frequency=0u"Hz", delay=0u"s", sys=Scanner())
    flip_angle, bandwidth, time_bw_product, phase_offset, freq_offset, delay = upreferred.((flip_angle, bandwidth, time_bw_product, phase_offset, freq_offset, delay)) .|> ustrip
    return PulseDesigner.block_pulse(FlipAngle(flip_angle), Bandwidth(bandwidth), TimeBwProduct(time_bw_product); phase_offset, freq_offset, delay, sys)
end

end
