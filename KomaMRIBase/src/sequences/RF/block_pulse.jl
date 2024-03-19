
function block_pulse(flip_angle::FlipAngle, duration::Duration;
    phase_offset=0, freq_offset=0, delay=0, sys=Scanner())

    flip_angle, duration = flip_angle.val, duration.val
    amplitude = flip_angle / (2π * γ * duration) * exp(im * phase_offset)
    dead_time = sys.RF_dead_time_T
    ring_down_time = sys.RF_ring_down_T
    delay = max(delay, dead_time)
    block_duration = delay + duration + ring_down_time

    rf = RF(amplitude, duration, freq_offset, delay)
    return Sequence([Grad(0, 0);;], [rf;;], [ADC(0, 0)], block_duration)
end

function block_pulse(flip_angle::FlipAngle, bandwidth::Bandwidth;
    phase_offset=0, freq_offset=0, delay=0, sys=Scanner())
    duration = 1 / (4 * bandwidth.val)
    return block_pulse(flip_angle, Duration(duration); phase_offset, freq_offset, delay, sys)
end

function block_pulse(flip_angle::FlipAngle, bandwidth::Bandwidth, time_bw_product::TimeBwProduct;
    phase_offset=0, freq_offset=0, delay=0, sys=Scanner())
    duration = time_bw_product.val / bandwidth.val
    return block_pulse(flip_angle, Duration(duration); phase_offset, freq_offset, delay, sys)
end

export block_pulse
