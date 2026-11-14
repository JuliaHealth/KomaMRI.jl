@testitem "Adaptive digital outputs and triggers" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    for types in ((1,), (2,), (1, 2))
        seq = Sequence()
        @addblock seq += Delay(1e-3)
        for type in types
            @addblock seq += Trigger(type, 1, 0.25e-3, 0.5e-3)
        end
        source = P.plot_seq(seq; adaptive=true).source
        visibility = [get(trace, :visible, true) for trace in source.blocks.data]
        embedded = P.plot_seq(seq)

        # Both event types retain absolute times, hover information and styling across zooms.
        expected_times = length(types) == 1 ? [1.25] : [1.25, 2.0]
        expected_labels = [type == 1 ? "Digital output · channel 1" : "Input trigger · channel 1" for type in types]
        @test source.triggers.times ≈ expected_times
        @test unique(last(embedded.data)[:text]) == expected_labels
        for interval in (source.blocks.full_range, [1.0, 1.5])
            payload = P.window_data(source, interval, 640, visibility)
            expected = filter(t -> interval[1] <= t <= interval[2], expected_times)
            @test [shape[:x0] for shape in payload.layout[:shapes]] ≈ expected
            @test unique(last(payload.data)[:x]) ≈ expected
            @test unique(last(payload.data)[:text]) == expected_labels[findall(t -> interval[1] <= t <= interval[2], expected_times)]
            @test last(payload.data)[:hovertemplate] == last(embedded.data)[:hovertemplate]
            @test first(payload.layout[:shapes])[:line] == first(embedded.layout[:shapes])[:line]
        end
    end
end
