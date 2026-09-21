@testitem "Slew-rate intervals" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    ramp = 1 / 128
    amplitudes = [1 / 64, -1 / 32, 1 / 128]
    triangle = Sequence(reshape([Grad(a, 0.0, ramp) for a in amplitudes], 3, 1))
    value_at(trace, time) = trace[:y][searchsortedfirst(trace[:x], time * 1e3)]
    source = plot_slew_rate(triangle; adaptive=true).source
    visibility = fill(true, length(source.data))
    detail = P.window_data(source, source.full_range, 640, visibility)

    # A triangle differentiates to constant positive/negative steps, including the first/last intervals.
    for traces in (plot_slew_rate(triangle).data, detail.data), axis in 1:3
        trace = traces[axis]
        @test trace[:mode] == "lines" && trace[:line][:shape] == "vh"
        @test first(trace[:x]) == 0 && last(trace[:x]) == 2ramp * 1e3
        @test value_at(trace, ramp / 2) == amplitudes[axis] / ramp
        @test value_at(trace, 3ramp / 2) == -amplitudes[axis] / ramp
    end

    pulses = 128
    pause = 1 / 2
    repeated = Sequence()
    for _ in 1:2
        for _ in 1:pulses
            @addblock repeated += triangle
        end
        @addblock repeated += Delay(pause)
    end
    source = plot_slew_rate(repeated; adaptive=true).source
    overview = P.window_data(source, source.full_range, 94, visibility)
    gap_start = pulses * 2ramp
    gap = [gap_start, gap_start + pause] .* 1e3
    detail = P.window_data(source, gap, 640, visibility)

    # Downsampling and zoom must keep zero slew throughout the delay between repetitions.
    @test overview.shown < overview.available
    for traces in (plot_slew_rate(repeated).data, overview.data, detail.data), axis in 1:3
        @test traces[axis][:line][:shape] == "vh"
        @test all(fraction -> value_at(traces[axis], gap_start + fraction * pause) == 0,
            (1 / 4, 1 / 2, 3 / 4))
    end
end
