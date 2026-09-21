@testitem "Adaptive phantom motion and spatial detail" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    grid = collect(range(-0.01, 0.01; length=41))
    x, y = repeat(grid, 41), repeat(grid; inner=41)
    density = collect(Float64, 1:length(x))
    obj = Phantom(; x, y, ρ=density,
        motion=translate(0.02, 0.0, 0.0, TimeRange(0.0, 1.0)))
    plot = P.plot_phantom_map(obj, :ρ; adaptive=true, time_samples=5)
    query = P.spatial_query(plot; width=40, height=40)
    # Camera arrays from browser messages are untyped, including on motion updates.
    query["projection"] = Any[query["projection"]...]
    overview = P.spatial_data(plot, query)

    # Reducing the screen cloud never changes its spin values, identity, or global color scale.
    trace = only(overview.data)
    ids = trace[:customdata]
    @test overview.shown < overview.available == length(obj)
    @test trace[:marker][:color] == density[ids]
    @test trace[:x] == 100 .* x[ids]
    @test trace[:marker][:cmax] == maximum(density)

    # Zooming retrieves original spins omitted in the overview, not interpolated replacements.
    query["width"] = query["height"] = 1000
    detail = P.spatial_data(plot, query)
    @test only(detail.data)[:customdata] == collect(eachindex(x))
    @test detail.shown == detail.available

    # Browser dimensions decoded as small unsigned integers must not overflow the pixel budget.
    query["width"], query["height"] = UInt16(575), UInt16(847)
    @test P.spatial_data(plot, query).shown == length(obj)

    # Time changes move the same IDs analytically; limits and the two-frame cache stay bounded.
    for frame in eachindex(plot.source.times)
        query["frame"] = frame
        moved = only(P.spatial_data(plot, query).data)
        moved_ids = moved[:customdata]
        expected = 100 .* (x[moved_ids] .+ 0.02 * plot.source.times[frame])
        @test moved[:x] == expected
        @test moved[:marker][:color] == density[moved_ids]
        @test moved[:marker][:cmax] == maximum(density)
        @test length(plot.source.frames) <= P.SPATIAL_FRAME_CACHE
    end

    # A spin that moves into a zoomed viewport is selected using its current, not reference, position.
    query["width"] = query["height"] = 40
    matrix = P.spatial_projection(plot)
    matrix[1, 1], matrix[1, 4] = 2 / 0.5, -2 * 2 / 0.5
    query["projection"] = vec(matrix)
    query["frame"] = 1
    @test P.spatial_data(plot, query).available == 0
    query["frame"] = length(plot.source.times)
    cropped = P.spatial_data(plot, query)
    visible = only(cropped.data)
    @test cropped.available > 0
    @test count(v -> 1.75 <= v <= 2.25, visible[:x]) == cropped.shown
    @test visible[:x] == 100 .* (x[visible[:customdata]] .+ 0.02)

    # A 3D translation preserves all three physical coordinates and IDs in a full-detail frame.
    volume = Phantom(x=[-0.01, 0.01], y=[0.01, -0.01], z=[0.0, 0.02],
        motion=translate(0.01, -0.02, 0.03, TimeRange(0.0, 1.0)))
    moving = P.plot_phantom_map(volume, :ρ; adaptive=true)
    query = P.spatial_query(moving)
    query["frame"] = length(moving.source.times)
    frame = only(P.spatial_data(moving, query).data)
    @test frame[:customdata] == [1, 2]
    for (axis, displacement) in zip((:x, :y, :z), (0.01, -0.02, 0.03))
        @test frame[axis] == 100 .* (getproperty(volume, axis) .+ displacement)
    end

    # A uniform 3D volume must retain samples throughout its depth, not only the first layer.
    grid = collect(-10:10) ./ 1000
    positions = vec(collect(Iterators.product(grid, grid, grid)))
    volume = Phantom(; x=first.(positions), y=getindex.(positions, 2), z=last.(positions))
    plot = P.plot_phantom_map(volume, :ρ; adaptive=true, view_2d=false)

    # Static 3D views do not reserve a bottom gutter; motion retains room for its time slider.
    @test plot.template.layout[:margin][:b] == 0
    @test moving.template.layout[:margin][:b] > 0

    query = P.spatial_query(plot; width=32, height=32)
    matrix = zeros(4, 4)
    matrix[1, 1] = matrix[2, 2] = matrix[4, 4] = 1
    query["projection"] = vec(matrix)
    for axis in 1:3
        permutation = circshift(collect(1:3), axis - 1)
        query["projection"] = vec(matrix[:, [permutation; 4]])
        reduced = P.spatial_data(plot, query)
        points = only(reduced.data)
        @test reduced.shown < reduced.available == length(volume)
        @test reduced.shown <= 32^2 ÷ P.SPATIAL_PIXEL_SIZE^2
        @test all((:x, :y, :z)) do dimension
            lower, upper = extrema(points[dimension])
            lower < -0.5 && upper > 0.5
        end
    end

    # Cropping removes outside spins, including when the visible cloud fits the pixel budget.
    matrix[1, 1] = matrix[2, 2] = 5
    query["projection"] = Any[vec(matrix)...]
    cropped = P.spatial_data(plot, query)
    points = only(cropped.data)
    @test all(x -> abs(x) <= 0.2, points[:x])
    @test count(zip(points[:x], points[:y])) do (x, y)
        abs(x) <= 0.2 && abs(y) <= 0.2
    end == cropped.shown
    @test points[:x] == 100 .* volume.x[points[:customdata]]
    query["width"] = query["height"] = 1000
    detail = P.spatial_data(plot, query)
    @test detail.shown == detail.available < length(volume)
    @test all(x -> abs(x) <= 0.2, only(detail.data)[:x])

    # Motion previews follow stable original spins; stopping restores omitted spins and their values.
    grid = collect(-12:12) ./ 1000
    positions = vec(collect(Iterators.product(grid, grid, grid)))
    volume = Phantom(; x=first.(positions), y=getindex.(positions, 2), z=last.(positions),
        ρ=collect(Float64, eachindex(positions)),
        motion=translate(0.01, 0.0, 0.0, TimeRange(0.0, 1.0)))
    moving = P.plot_phantom_map(volume, :ρ; adaptive=true, view_2d=false)
    query = P.spatial_query(moving)
    query["preview"] = true
    preview = P.spatial_data(moving, query)
    ids = only(preview.data)[:customdata]
    @test preview.shown == length(ids) < preview.available == length(volume)
    @test length(ids) <= P.SPATIAL_MOTION_BUDGET
    for frame in eachindex(moving.source.times)
        query["frame"] = frame
        preview_trace = only(P.spatial_data(moving, query).data)
        @test preview_trace[:customdata] == ids
        @test preview_trace[:x] == 100 .* (volume.x[ids] .+ 0.01 * moving.source.times[frame])
        @test preview_trace[:marker][:color] == volume.ρ[ids]
    end
    query["preview"] = false
    detail = P.spatial_data(moving, query)
    @test detail.complete && detail.shown == detail.available == length(volume)
    @test only(detail.data)[:customdata] == collect(eachindex(volume.x))

    # Preview budgets refill from the zoomed region, revealing spins omitted at full FOV.
    query["preview"] = true
    matrix = zeros(4, 4)
    matrix[1, 1], matrix[2, 2], matrix[4, 4] = 5, 5, 1
    matrix[1, 4] = -5
    query["projection"] = vec(matrix)
    zoomed = P.spatial_data(moving, query)
    zoom_ids = only(zoomed.data)[:customdata]
    @test zoomed.shown == zoomed.available < P.SPATIAL_MOTION_BUDGET
    @test !isempty(setdiff(zoom_ids, ids))
    @test all(x -> 0.8 <= x <= 1.2, only(zoomed.data)[:x])
end

@testitem "Phantom speed and reinjection" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    obj = Phantom(x=[0.0, 0.0], motion=flowpath(
        [0.0 0.03 0.06; 0.0 1.0 1.03], [0.0 0.04 0.08; 0.0 1.0 1.04],
        zeros(2, 3), [false false false; false true false], TimeRange(0.0, 2.0)))
    plot = P.plot_phantom_map(obj, :speed; adaptive=true, time_samples=3)
    query = P.spatial_query(plot)

    # Exit is shown at the pre-reset position; entry is shown at the new position one frame later.
    speed, exiting, entering = P.spatial_data(plot, query).data
    @test exiting[:customdata] == [2]
    @test (exiting[:x], exiting[:y]) == ([0.0], [0.0])
    @test exiting[:marker][:color] == "#ff40df"
    @test isempty(entering[:customdata])
    query["frame"] = 2
    speed, exiting, entering = P.spatial_data(plot, query).data
    @test isempty(exiting[:customdata])
    @test entering[:customdata] == [2]
    @test (entering[:x], entering[:y]) == ([100.0], [100.0])
    @test entering[:marker][:color] == "#00e5ff"

    # A 3-4-5 cm displacement has speed 5 cm/s; reinjection is not a high velocity.
    @test speed[:marker][:color] == [5.0]
    @test speed[:marker][:cmax] == 5.0
    query["frame"] = 3
    speed, exiting, entering = P.spatial_data(plot, query).data
    @test speed[:marker][:color] == [5.0, 5.0]
    @test isempty(exiting[:customdata]) && isempty(entering[:customdata])

    # A periodic file's final frame shows the recorded endpoint, not an invented all-spin reset.
    obj.motion.time = TimeCurve(t=[0.0, 1.0, 2.0], t_unit=[0.0, 0.5, 1.0], periodic=true)
    looping = P.plot_phantom_map(obj, :speed; adaptive=true, time_samples=3)
    @test first(looping.template.data)[:marker][:cmax] == 5.0
    last_frame = P.spatial_data(looping, query).data
    @test last_frame[1][:x] == [6.0, 103.0]
    @test last_frame[1][:marker][:color] == [5.0, 5.0]
    @test isempty(last_frame[2][:customdata]) && isempty(last_frame[3][:customdata])
    @test all(isnothing, last_frame[2][:x]) && all(isnothing, last_frame[3][:x])
    @test obj.motion.time.periodic
    obj.motion.time = TimeRange(0.0, 2.0)

    # Portable frames retain velocities and both event overlays, switching all three traces together.
    portable = P.plot_phantom_map(obj, :speed; time_samples=3)
    @test portable.data[2][:customdata] == [2]
    @test portable.data[4][:marker][:color] == [5.0]
    @test portable.data[6][:customdata] == [2]
    @test findall(portable.layout[:sliders][1][:steps][2][:args][1][:visible]) == 4:6
end

@testitem "Phantom path frame interpolation" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    # Frame evaluation matches simulation interpolation at knots, between knots, and outside the motion.
    for T in (Float32, Float64), spins in (AllSpins(), SpinRange(1:2:3), SpinRange(2:2))
        selected = collect(1:3)[KomaMRIBase.get_indexing_range(spins)]
        displacement = T[0 0.125 0.5; 0 -0.25 -0.5; 0 0.5 0.25][selected, :]
        obj = Phantom(x=T[0, 0.125, 0.25], motion=Motion(
            Path(displacement, -displacement, 2displacement), TimeRange(T(1), T(3)), spins))
        for t in T[0, 1, 1.375, 2, 2.25, 3, 4]
            @test P.phantom_coordinates(obj.motion, obj, t) ==
                get_spin_coords(obj.motion, obj.x, obj.y, obj.z, [t]')
        end
    end
end

@testitem "Phantom velocity components" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    step = 1 / 32
    displacement = step .* [0 1 3; 0 -1 -3; 0 0 0; 0 0 0]
    displacement[4, 2:3] .= 1
    resets = falses(4, 3)
    resets[4, 2] = true
    obj = Phantom(x=zeros(4), motion=flowpath(displacement, 2displacement,
        -4displacement, resets, TimeRange(0.0, 2.0)))

    for (key, velocity) in zip((:vx, :vy, :vz), 100step .* (1, 2, -4))
        plot = P.plot_phantom_map(obj, key; adaptive=true, time_samples=3)
        query = P.spatial_query(plot)

        # Forward differences at the first frame and backward differences agree in cm/s.
        for frame in 1:2
            query["frame"] = frame
            moving, exiting, entering = P.spatial_data(plot, query).data
            @test moving[:marker][:color] == [velocity, -velocity, 0]
            @test (frame == 1 ? exiting : entering)[:customdata] == [4]
        end

        # The scale spans the whole motion symmetrically; reinjection jumps do not enlarge it.
        marker = first(plot.template.data)[:marker]
        @test (marker[:cmin], marker[:cmax]) == (-2abs(velocity), 2abs(velocity))
        @test marker[:colorscale] == [[0, "blue"], [0.5, "white"], [1, "red"]]

        # Later motion changes magnitude without rescaling, including in the portable plot.
        query["frame"] = 3
        moving, exiting, entering = P.spatial_data(plot, query).data
        @test moving[:marker][:color] == [2velocity, -2velocity, 0, 0]
        @test isempty(exiting[:customdata]) && isempty(entering[:customdata])
        portable = P.plot_phantom_map(obj, key; time_samples=3)
        @test portable.data[7][:marker][:color] == [2velocity, -2velocity, 0, 0]
    end

    # A reduced flow preview keeps original spin IDs, analytical velocities, and reset identities.
    count = 2P.SPATIAL_MOTION_BUDGET
    step = collect(1:count) ./ 2^20
    offsets = hcat(zeros(count), step, 2step)
    resets = falses(count, 3)
    resets[2:2:end, 2] .= true
    obj = Phantom(x=zeros(count), motion=flowpath(offsets, zeros(count, 3), zeros(count, 3),
        resets, TimeRange(0.0, 2.0)))
    plot = P.plot_phantom_map(obj, :vx; adaptive=true, time_samples=3)
    query = P.spatial_query(plot)
    query["preview"], query["frame"] = true, 2
    moving, exiting, entering = P.spatial_data(plot, query).data
    @test length([moving[:customdata]; entering[:customdata]]) == P.SPATIAL_MOTION_BUDGET < count
    @test all(isodd, moving[:customdata]) && all(iseven, entering[:customdata])
    @test moving[:marker][:color] == 100step[moving[:customdata]]
    @test entering[:x] == 100step[entering[:customdata]]

    # Strided motion groups remain correct even when random preview IDs cannot form a SpinRange.
    obj.motion = Motion(Translate(0.01, 0.0, 0.0), TimeRange(0.0, 1.0), SpinRange(1:2:count))
    plot = P.plot_phantom_map(obj, :ρ; adaptive=true, time_samples=2)
    query = P.spatial_query(plot)
    query["preview"], query["frame"] = true, 2
    trace = only(P.spatial_data(plot, query).data)
    @test trace[:x] == Float64.(isodd.(trace[:customdata]))
end

@testitem "Shared phantom properties" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P
    obj = Phantom(x=[0.0, 0.01], ρ=[1.0, 2.0], T1=[1.0, 2.0],
        motion=translate(0.02, -0.01, 0.0, TimeRange(0.0, 1.0)))
    plot = P.plot_phantom(obj; key=:vx, adaptive=true, time_samples=3)
    query = P.spatial_query(plot)
    query["frame"] = 2
    velocity = P.spatial_data(plot, query)
    coordinates = P.phantom_frame(plot.source, 2)

    # Property switches keep the selected time, original spin IDs and cached geometry.
    query["property"] = "T1"
    relaxation = P.spatial_data(plot, query)
    @test first(velocity.data)[:customdata] == first(relaxation.data)[:customdata] == [1, 2]
    @test first(velocity.data)[:x] == first(relaxation.data)[:x]
    @test P.phantom_frame(plot.source, 2) === coordinates
    @test first(velocity.data)[:marker][:color] == [2.0, 2.0]
    @test first(relaxation.data)[:marker][:color] == [1000.0, 2000.0]
    @test first(relaxation.data)[:marker][:colorbar][:ticksuffix] == " ms"

    # Portable property buttons recolor existing time traces without replacing coordinates or visibility.
    portable = P.plot_phantom(obj; time_samples=3)
    buttons = [button for menu in portable.layout[:updatemenus] for button in menu[:buttons]]
    # One row separates intrinsic properties from motion-only velocity controls.
    @test [b[:label] for b in buttons] == ["ρ", "T1", "T2", "T2s", "Δw", "|v|", "vx", "vy", "vz"]
    @test all(menu -> menu[:y] == 1 && menu[:pad][:b] > portable.layout[:title][:pad][:b], portable.layout[:updatemenus])
    @test only(portable.layout[:annotations])[:text] == "—"
    vx = only(b for b in buttons if b[:label] == "vx")[:args][1]
    @test !any(haskey(vx, key) for key in ("x", "y", "z", "visible"))
    @test vx["marker.color"][1] == vx["marker.color"][4] == [2.0, 2.0]
    @test findall(portable.layout[:sliders][1][:steps][2][:args][1][:visible]) == 4:6

    # A stationary phantom still supports property controls without inventing a motion slider.
    stationary = P.plot_phantom(Phantom(x=[0.0, 0.01], T1=[1.0, 2.0]))
    @test !only(stationary.layout[:sliders])[:visible]
    @test [b[:label] for b in only(stationary.layout[:updatemenus])[:buttons]] == ["ρ", "T1", "T2", "T2s", "Δw"]
    @test isempty(stationary.layout[:annotations])
    t1 = only(b for menu in stationary.layout[:updatemenus] for b in menu[:buttons] if b[:label] == "T1")
    @test t1[:args][1]["marker.color"][1] == [1000.0, 2000.0]
end

@testitem "Adaptive coil grid and components" tags=[:plots] begin
    using KomaMRIBase
    import KomaMRIPlots as P

    coordinates = [-0.02, 0.0, 0.02]
    receiver = ArbitraryCoilSens(coordinates, coordinates, [0.0],
        cat(fill(2.0im, 3, 3, 1, 1), fill(-4.0 + 0im, 3, 3, 1, 1); dims=4))
    sys = Scanner(; receiver)
    plot = P.plot_coil_sens(sys; adaptive=true, spacing=0.001)
    query = P.spatial_query(plot; width=40, height=40)
    query["projection"] = Any[query["projection"]...]
    overview = P.spatial_data(plot, query)

    # Coarse samples remain on the requested physical grid and retain the true complex gain.
    @test overview.coarse
    @test overview.shown < prod(length, plot.source.grid)
    query["width"] = query["height"] = 1000
    full = P.spatial_data(plot, query)
    @test !full.coarse
    @test full.shown == prod(length, plot.source.grid)
    @test sort(unique(only(full.data)[:x])) == 100 .* collect(plot.source.grid[1])
    # Linear interpolation of a constant field can round its weighted terms in Float64.
    @test maximum(abs, only(full.data)[:marker][:color] .- 2) <= 8eps(2.0)
    @test full.coloraxis[:cmax] == 4

    # Coil and phase selection commute with spatial refinement; phase has a fixed cyclic scale.
    query["coil"], query["component"] = 2, 2
    phase = P.spatial_data(plot, query)
    @test all(==(Float64(π)), only(phase.data)[:marker][:color])
    @test (phase.coloraxis[:cmin], phase.coloraxis[:cmax]) == (-π, π)
    query["component"] = 1
    @test P.spatial_data(plot, query).coloraxis[:cmax] == 4

    # A fine 3D grid is implicit: displayed work is bounded, and default plots remain standalone.
    volume = P.plot_coil_sens(Scanner(); adaptive=true, spacing=0.0001)
    reduced = P.spatial_data(volume, P.spatial_query(volume; width=100, height=100))
    @test prod(length, volume.source.grid) > 1_000_000
    @test reduced.coarse && reduced.shown <= 100 * 100 ÷ P.SPATIAL_PIXEL_SIZE^2
    @test all(==(1), only(reduced.data)[:marker][:color])
    @test P.plot_coil_sens(sys) isa P.PlotlyBase.Plot
    @test P.plot_phantom_map(Phantom(x=[0.0]), :ρ) isa P.PlotlyBase.Plot
end
