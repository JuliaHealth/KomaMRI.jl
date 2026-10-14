@testitem "Phantom slabs follow current spin positions" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    spacing = 1 / 128
    obj = Phantom(x=[-spacing, 0, spacing], y=[-spacing, 0, spacing],
        z=[0, spacing, 2spacing], ρ=[1.0, 2.0, 3.0],
        motion=translate(spacing, 0.0, 0.0, TimeRange(0.0, 1.0)))
    live = P.plot_phantom_slice(obj; thickness=spacing, center=(0, 0, 0),
        adaptive=true, time_samples=2)
    query = P.spatial_query(live)
    selected(data, panel) = reduce(vcat,
        (trace[:customdata] for trace in data if
            startswith(get(trace, :uid, ""), "phantom-slice-$panel-")); init=Int[])

    # Fixed slabs select moving spins by their current position, without losing identity or density.
    for (frame, expected) in ((1, [[1], [2], [2]]), (2, [[1], [1], [2]]))
        query["frame"] = frame
        result = P.spatial_data(live, query)
        @test [selected(result.data, panel) for panel in 1:3] == expected
        @test [count.shown for count in result.slice_counts] == length.(expected)
        for trace in result.data[1:end-1]
            @test trace[:marker][:color] == obj.ρ[trace[:customdata]]
        end
    end

    # Offline slices retain the same initial selection and the complete analytical motion frames.
    portable = P.plot_phantom_slice(obj; thickness=spacing, center=(0, 0, 0), time_samples=2)
    @test [selected(portable.data, panel) for panel in 1:3] == [[1], [2], [2]]
    frames = portable.layout[:meta]["phantom_slices"].frames
    @test [frame.coordinates[1] for frame in frames] ==
        [100 .* obj.x, 100 .* (obj.x .+ spacing)]

    # Empty slabs retain the phantom's color scale rather than dropping the colorbar.
    query["slices"]["center"] = fill(100.0, 3)
    empty_slabs = P.spatial_data(live, query)
    @test all(count -> count.shown == count.available == 0, empty_slabs.slice_counts)
    colorbar = only(filter(trace -> trace[:uid] == "phantom-slice-colorbar", empty_slabs.data))
    @test colorbar[:marker][:showscale] && colorbar[:marker][:color] == [0, maximum(obj.ρ)]
end
