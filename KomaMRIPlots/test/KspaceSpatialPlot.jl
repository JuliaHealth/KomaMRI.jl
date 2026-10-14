@testitem "Adaptive k-space segments and ADC provenance" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    seq = PulseDesigner.EPI_example()
    # A single-point origin must not create a zero-length GL line with an invalid hover pick.
    @test P.plot_kspace(seq).data[3][:mode] == "markers"
    @test P.plot_kspace(seq; adaptive=true).template.data[3][:mode] == "markers"
    plot = P.plot_kspace(seq; adaptive=true, view_2d=true)
    query = P.spatial_query(plot; width=20, height=20)
    payload = P.spatial_data(plot, query)
    line, adc = payload.data[1:2]

    # Decimation retains measured coordinates/times and cannot join omitted trajectory pieces.
    ids = adc[:customdata]
    @test payload.shown < payload.available
    @test adc[:x] == plot.source.adc[ids, 1]
    @test adc[:text] == round.(plot.source.times[ids] .* 1e3; digits=3)
    @test adc[:marker][:color] == ids .- 1
    @test length(ids) <= P.KSPACE_DISPLAY_BUDGET
    for (j, i) in enumerate(line[:customdata])
        @test line[:x][j] == (iszero(i) ? nothing : plot.source.trajectory[i, 1])
        j == 1 || iszero(i) || iszero(line[:customdata][j-1]) ||
            @test i == line[:customdata][j-1] + 1
    end
    @test !isempty(P.PlotlyBase.JSON.json(payload))

    # A fully visible trajectory keeps every original vertex without duplicated endpoints.
    full = P.spatial_data(plot, P.spatial_query(plot; width=900, height=900))
    @test full.data[1][:customdata] == collect(1:size(plot.source.trajectory, 1))
    @test full.data[1][:x] == plot.source.trajectory[:, 1]

    # Hiding either layer omits its data and counts, without changing the visible layer.
    query["visible"] = [false, true, true]
    hidden = P.spatial_data(plot, query)
    @test isempty(hidden.data[1][:x])
    @test hidden.data[2][:customdata] == adc[:customdata]
    @test hidden.shown == length(adc[:customdata])
    query["visible"] = [true, false, true]
    hidden = P.spatial_data(plot, query)
    @test isempty(hidden.data[2][:x])
    @test hidden.data[1][:customdata] == line[:customdata]
    delete!(query, "visible")

    # A segment crossing a zoomed view remains visible even when both endpoints are outside.
    matrix = Float64[i == j for i in 1:4, j in 1:4]
    crossing = [-2.0 0 0; 2 0 0; 2 2 0]
    segments, available = P.kspace_segments(crossing, matrix, 100, 100, 100)
    @test segments == [1] && available == 1

    # Zooming into an individual ADC restores its exact sample and excludes distant ADCs.
    target = size(plot.source.adc, 1) ÷ 2
    matrix[1, 1] = matrix[2, 2] = 100
    matrix[1:2, 4] .= -100 .* plot.source.adc[target, 1:2]
    matrix[3, :] .= 0
    query["projection"] = vec(matrix)
    query["width"] = query["height"] = 900
    zoom = P.spatial_data(plot, query)
    ids = zoom.data[2][:customdata]
    # A sample keeps its acquisition-order color and global scale when the viewport changes.
    @test zoom.data[2][:marker][:color] == ids .- 1
    for view_2d in (false, true)
        marker = P.plot_kspace(seq; view_2d).data[2][:marker]
        @test marker[:color] == collect(0:length(plot.source.times)-1)
        for key in (:colorscale, :cmin, :cmax, :showscale, :colorbar)
            @test marker[key] == adc[:marker][key] == zoom.data[2][:marker][key]
        end
    end
    # The refined viewport restores the target sample without retaining distant samples.
    @test target in ids
    @test all(i -> all(abs.(plot.source.adc[i, 1:2] - plot.source.adc[target, 1:2]) .<= 0.01), ids)
    @test length(ids) < size(plot.source.adc, 1)
end
