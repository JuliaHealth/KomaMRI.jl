
function block_pulse_duration(x::Unitful.Quantity)
    println("Hola")
end

function block_pulse(
    flip_angle::Angle, duration::Unitful.Time;
    phase_offset::Angle=0u"°", freq_offset::Frequency=0u"Hz", delay=0u"s", sys=Scanner()
    )

    # Compute amplitude
    flip_angle_rad = flip_angle |> u"rad" |> ustrip
    amplitude = flip_angle_rad / (2π * uγ * duration) * exp(im * phase_offset)

    # Give units to hardware constrainsts
    dead_time = sys.RF_dead_time_T * u"s"
    ring_down_time = sys.RF_ring_down_T * u"s"

    # Compute delay and block duration
    delay = max(delay, dead_time)
    dur = delay + duration + ring_down_time

    # Get values in SI units
    amplitude, duration, freq_offset, delay, dur = upreferred.((amplitude, duration, freq_offset, delay, dur))

    # Return the sequence
    rf = RF(amplitude, duration, freq_offset, delay)
    return Sequence([Grad(0, 0);;], [rf;;], [ADC(0, 0)], dur)

end

function block_pulse(
    flip_angle::Angle, bandwidth::Unitful.Frequency;
    phase_offset::Angle=0u"°", freq_offset::Frequency=0u"Hz", delay=0u"s", sys=Scanner()
    )

    println("Frequency")
    duration = 1 / (4 * bandwidth)
    amplitude = flip_angle / (2π * γ * duration)
    # ...

end
function block_pulse(
    flip_angle::Angle, bandwidth::DimensionlessQuantity;
    phase_offset::Angle=0u"°", freq_offset::Frequency=0u"Hz", delay=0u"s", sys=Scanner()
    )

    println("TBP")
    duration = 1 / (4 * bandwidth)
    amplitude = flip_angle / (2π * γ * duration)
    # ...
end

"""
"""
function block_pulse(flip_angle; duration=0, freq_offset=0, phase_offset=0,
    time_bw_product=0, bandwidth=0, delay=0,
    dead_time=0, ring_down_time=0, sys=nothing)

    if !isnothing(sys)
        dead_time = sys.RF_dead_time_T
        ring_down_time = sys.RF_ring_down_T
    end

    if duration == 0
        if time_bw_product > 0
            duration = time_bw_product / bandwidth
        elseif bandwidth > 0
            duration = 1 / (4 * bandwidth)
        else
            @error "Either bandwidth or duration must be defined"
        end
    end
    amplitude = flip_angle / (2π * γ * duration) * cis(phase_offset)

    if dead_time > delay
        delay = dead_time
    end

    dly = Delay(0)
    if ring_down_time > 0
        dly = Delay(delay + duration + ring_down_time)     # I NEED a review
    end

    return RF(amplitude, duration, freq_offset, delay), dly
end
