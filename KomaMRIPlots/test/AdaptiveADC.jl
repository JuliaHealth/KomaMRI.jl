@testitem "Adaptive ADC windows and samples" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    seq = Sequence()
    for _ in 1:128
        @addblock seq += ADC(16, 15e-6, 10e-6)
        @addblock seq += Delay(1e-3)
    end
    source = P.plot_seq(seq; adaptive=true).source
    visibility = [get(trace, :visible, true) for trace in source.blocks.data]
    overview = P.window_data(source, source.blocks.full_range, 142, visibility)
    window, markers = overview.data[8:9]

    # Retained acquisition windows remain rectangular; dense ADC markers are hidden, not subsampled.
    x, y = window[:x], window[:y]
    @test all(2:length(x)) do i
        isnothing(x[i-1]) || isnothing(x[i]) || x[i-1] == x[i] || y[i-1] == y[i]
    end
    @test isempty(markers[:x])

    # At about 2.3 px spacing, every sample appears at the simulator's original acquisition time.
    interval = [0.0, 0.25]
    detail = P.window_data(source, interval, 640, visibility)
    markers = detail.data[9]
    times = filter(t -> !isnothing(t) && interval[1] <= t <= interval[2], markers[:x])
    @test times == get_adc_sampling_times(seq)[1:16] .* 1e3
    @test all(==(1), filter(!isnothing, markers[:y]))
    @test detail.shown == detail.available

    # Default ADC samples retain the same vertical tick at different zoom levels.
    @test markers[:marker][:symbol] == "line-ns"
    closer = P.window_data(source, [0.0, 0.05], 640, visibility)
    @test closer.data[9][:marker][:symbol] == "line-ns"

    # Single-sample acquisitions mark their centre and obey the same zoom-density rule.
    single = Sequence()
    for _ in 1:128
        @addblock single += ADC(1, 10e-6)
    end
    source = P.plot_seq(single; adaptive=true).source
    overview = P.window_data(source, source.blocks.full_range, 142, visibility)
    @test isempty(overview.data[9][:x])
    detail = P.window_data(source, [0.0, dur(first(single.ADC)) * 1e3], 640, visibility)
    @test first(detail.data[9][:x]) == only(get_adc_sampling_times(single[1:1])) * 1e3
end

@testitem "Plot keyword modes" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    seq = Sequence()
    @addblock seq += ADC(16, 15e-6)

    # Nonadaptive-only keywords cannot silently switch live sequences to full trace construction.
    for options in ((; gl=true), (; show_adc=true), (; max_rf_samples=10),
        (; show_seq_blocks=true), (; freq_in_phase=true), (; show_rf_frame=true),
        (; xaxis="x2"), (; yaxis="y2"), (; showlegend=false))
        @test_throws ArgumentError P.plot_seq(seq; adaptive=true, options...)
    end
    portable = P.plot_seq(seq; show_adc=true)
    adc = only(filter(trace -> trace[:name] == "ADC", portable.data))
    @test collect(skipmissing(adc[:x])) == get_adc_sampling_times(seq) .* 1e3

    # Fixed spin caps remain effective for embedded plots and are rejected for both live phantom APIs.
    obj = Phantom(x=collect(range(-0.01, 0.01; length=21)))
    for plotter in (obj -> P.plot_phantom(obj; max_spins=5),
        obj -> P.plot_phantom_map(obj, :ρ; max_spins=5))
        @test length(first(plotter(obj).data)[:x]) <= 5
    end
    @test_throws ArgumentError P.plot_phantom(obj; adaptive=true, max_spins=5)
    @test_throws ArgumentError P.plot_phantom_map(obj, :ρ; adaptive=true, max_spins=5)
end
