"""
Returns a boolean value indicating whether two vectors are approximately equal.
"""
function safe_approx(v1, v2)
    if length(v1) != length(v2)
        return false
    end
    return v1 ≈ v2
end

"""
Returns a boolean value indicating whether two vector of angles are equal.
"""
function safe_approx_angles(a1, a2)
    if length(a1) != length(a2)
        return false
    end
    return abs(sum(exp.(1im * 2π * (a1)) .* exp.(-1im * 2π * (a2)))) / length(a1) ≈ 1
end

"""
Returns a boolean value indicating whether a vector is present in a list of vectors.
"""
function not_in_list(vec, vec_list)
    for arr in vec_list
        if safe_approx(vec, arr)
            return false
        end
    end
    return true
end

"""
Returns a boolean value indicating whether a event is present in a list of events.
"""
function not_in_list_event(event, event_list)
    for obj in event_list
        if event ≈ obj
            return false
        end
    end
    return true
end

"""
Returns a boolean value indicating whether a vector is present in a list of vectors.
"""
function not_in_list_angles(angle, angle_list)
    for phase in angle_list
        if safe_approx_angles(angle, phase)
            return false
        end
    end
    return true
end

"""
Returns the events and IDs which are unique given an input vector of events.
"""
function get_unique_events(events::Vector)
    events_obj, events_id, id_cnt = [], [], 1
    for obj in events
        if not_in_list_event(obj, events_obj)
            push!(events_obj, obj)
            push!(events_id, id_cnt)
            id_cnt += 1
        end
    end
    return (obj = events_obj , id = events_id)
end

"""
Returns the vector of IDs for all the blocks of an array of events.
It is neccessary to add the input vector `events` which contains the uniques objects
(RF, Grad, or ADC) with its repective ID.
"""
function get_events_id(event_array, events)
    events_id = fill(0, length(event_array))
    for (blk, eve) in enumerate(event_array)
        for (obj, id) in zip(events.obj, events.id)
            if eve ≈ obj
                events_id[blk] = id
                break
            end
        end
    end
    return events_id
end

"""
Separates the `grs` into the gradients with arbitrary waveform and trapezoidal
"""
function get_unique_gradients(grs)
    grads = (obj = [], id = [])
    traps = (obj = [], id = [])
    for (obj, id) in zip(grs.obj, grs.id)
        if length(obj.A) == 1
            push!(traps.obj, obj)
            push!(traps.id, id)
        else
            push!(grads.obj, obj)
            push!(grads.id, id)
        end
    end
    return grads, traps
end

"""
Returns the unique shapes for the magnitude, angle and time of the `rfs` vector.
Requires an initial integer counter `id_shape_cnt` to asign IDs incrementally.
"""
function get_rf_shapes(rfs, id_shape_cnt, Δt_rf)
    # Find the unique shapes (magnitude, phase and time shapes) and assign IDs
    rfs_abs, rfs_ang, rfs_tim = [], [], []
    rfs_abs_id, rfs_ang_id, rfs_tim_id = [], [], []
    for obj in rfs.obj
        shape_abs = abs.(obj.A) / maximum(abs.(obj.A))
        if not_in_list(shape_abs, rfs_abs)
            push!(rfs_abs, shape_abs)
            push!(rfs_abs_id, id_shape_cnt)
            id_shape_cnt += 1
        end
        shape_ang = mod.(angle.(obj.A), 2π) / 2π
        ang = shape_ang .- shape_ang[1]
        list_ang_unique = [shape .- shape[1] for shape in rfs_ang]
        if not_in_list_angles(ang, list_ang_unique)
            push!(rfs_ang, shape_ang)
            push!(rfs_ang_id, id_shape_cnt)
            id_shape_cnt += 1
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / Δt_rf
            if not_in_list(shape_tim, rfs_tim)
                push!(rfs_tim, shape_tim)
                push!(rfs_tim_id, id_shape_cnt)
                id_shape_cnt += 1
            end
        end
    end
    rfmags = (data = rfs_abs, id = rfs_abs_id)
    rfangs = (data = rfs_ang, id = rfs_ang_id)
    rftimes = (data = rfs_tim, id = rfs_tim_id)
    return rfmags, rfangs, rftimes, id_shape_cnt
end

"""
Returns the unique shapes for the amplitude and time of the `grads` vector.
Requires an initial integer counter `id_shape_cnt` to asign IDs incrementally.
"""
function get_grad_shapes(grads, id_shape_cnt, Δt_gr)
    # Find shapes for magnitude and time gradients
    grads_amp, grads_tim = [], []
    grads_amp_id, grads_tim_id = [], []
    for obj in grads.obj
        shape_amp = obj.A / maximum(abs.(obj.A))
        if not_in_list(shape_amp, grads_amp)
            push!(grads_amp, shape_amp)
            push!(grads_amp_id, id_shape_cnt)
            id_shape_cnt = id_shape_cnt + 1
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / Δt_gr
            if not_in_list(shape_tim, grads_tim)
                push!(grads_tim, shape_tim)
                push!(grads_tim_id, id_shape_cnt)
                id_shape_cnt = id_shape_cnt + 1
            end
        end
    end
    gramps = (data = grads_amp, id = grads_amp_id)
    grtimes = (data = grads_tim, id = grads_tim_id)
    return gramps, grtimes, id_shape_cnt
end

"""
Defines the library to be written in the [BLOCKS] section
Columns of block_events: [blk, id_rf, id_gx, id_gy, id_gz, id_adc, id_ext]
"""
function get_block_events(seq, rfs, grs, adcs)
    r = get_events_id(seq.RF, rfs)
    x = get_events_id(seq.GR.x, grs)
    y = get_events_id(seq.GR.y, grs)
    z = get_events_id(seq.GR.z, grs)
    a = get_events_id(seq.ADC, adcs)
    block_events = [[b, 0, r[b], x[b], y[b], z[b], a[b], 0] for (b, _) in enumerate(seq)]
    for row in block_events
        blk = row[1]
        bd = seq.DUR[blk] / seq.DEF["BlockDurationRaster"]
        bdr = round(bd)
        if abs(bdr - bd) > 1e-6
            @warn "Block $blk duration rounded"
        end
        row[2] = bdr
    end
    return block_events
end

"""
Defines the library to be written in the [RF] section
Columns of rf_library: [id, amp, id_mag, id_phase, id_time, delay, freq, phase]
"""
function get_rf_library(rfs, rfmags, rfangs, rftimes, Δt_rf)
    rf_library = [[id, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for id in rfs.id]
    for (i, row) in enumerate(rf_library)
        obj = rfs.obj[i]
        row[2] = γ * maximum(abs.(obj.A))
        row[7] = obj.Δf
        shape_abs = abs.(obj.A) / maximum(abs.(obj.A))
        for (shape_abs_unique, id_abs) in zip(rfmags.data, rfmags.id)
            if safe_approx(shape_abs, shape_abs_unique)
                row[3] = id_abs
                break
            end
        end
        shape_ang = mod.(angle.(obj.A), 2π) / 2π
        ang = shape_ang .- shape_ang[1]
        for (shape_ang_unique, id_ang) in zip(rfangs.data, rfangs.id)
            ang_unique = shape_ang_unique .- shape_ang_unique[1]
            if safe_approx_angles(ang, ang_unique)
                row[4] = id_ang
                row[8] = angle(
                    sum(exp.(1im * 2π * shape_ang) .* exp.(-1im * 2π * shape_ang_unique)) /
                    length(shape_ang),
                )
                break
            end
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / Δt_rf
            for (shape_tim_unique, id_tim) in zip(rftimes.data, rftimes.id)
                if safe_approx(shape_tim, shape_tim_unique)
                    row[5] = id_tim
                    break
                end
            end
        end
        delay_compensation_rf_koma = (row[5] == 0) * Δt_rf / 2
        row[6] = round((obj.delay - delay_compensation_rf_koma) / Δt_rf) * Δt_rf * 1e6
    end
    return rf_library
end

"""
Defines the library to be written in the [GRADIENTS] section
Columns of grad_library_arb: [id, amp, id_amp, id_time, delay]
"""
function get_grad_library_arb(grads, gramps, grtimes, Δt_gr)
    grad_library_arb = [[id, 0.0, 0.0, 0.0, 0.0] for id in grads.id]
    for (i, row) in enumerate(grad_library_arb)
        obj = grads.obj[i]
        row[2] = γ * maximum(abs.(obj.A))    # this always stores positive values, the waveform vector have the respective positive or negative values
        row[5] = round(1e6 * obj.delay)
        shape_amp = obj.A / maximum(abs.(obj.A))
        for (shape_amp_unique, id_amp) in zip(gramps.data, gramps.id)
            if safe_approx(shape_amp, shape_amp_unique)
                row[3] = id_amp
                break
            end
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / Δt_gr
            for (shape_tim_unique, id_tim) in zip(grtimes.data, grtimes.id)
                if safe_approx(shape_tim, shape_tim_unique)
                    row[4] = id_tim
                    break
                end
            end
        end
    end
    return grad_library_arb
end

"""
Defines the library to be written in the [TRAP] section
Columns of grad_library_trap: [id, amp, rise, flat, fall, delay]
"""
function get_grad_library_trap(traps)
    grad_library_trap = [[id, 0.0, 0.0, 0.0, 0.0, 0.0] for id in traps.id]
    for (i, row) in enumerate(grad_library_trap)
        obj = traps.obj[i]
        row[2] = γ * obj.A
        row[3] = 1e6 * obj.rise
        row[4] = 1e6 * obj.T
        row[5] = 1e6 * obj.fall
        row[6] = 1e6 * obj.delay
    end
    return grad_library_trap
end

"""
Defines the library to be written in the [ADC] section
Columns of adc_library: [id, num, dwell, delay, freq, phase]
"""
function get_adc_library(adcs)
    adc_library = [[id, 0.0, 0.0, 0.0, 0.0, 0.0] for id in adcs.id]
    for (i, row) in enumerate(adc_library)
        obj = adcs.obj[i]
        row[2] = obj.N
        row[3] = obj.T * 1e9 / (obj.N - 1)
        row[4] = (obj.delay - 0.5 * obj.T / (obj.N - 1)) * 1e6
        row[5] = obj.Δf
        row[6] = obj.ϕ
    end
    return adc_library
end

"""
Defines the library to be written for the [SHAPES] section
Elements of shape_library: [id, num, data]
"""
function get_shape_library(rfmags, rfangs, rftimes, gramps, grtimes)
    shape_library = []
    for shapes in [rfmags, rfangs, rftimes, gramps, grtimes]
        for (data, id) in zip(shapes.data, shapes.id)
            num, data = compress_shape(data)
            push!(shape_library, (id = id, num = num, data = data))
        end
    end
    return shape_library
end

"""
Returns the with its libraries to be written in the obj file
"""
function get_pulseq_object(seq::Sequence)
    # Get the raster times
    Δt_rf = seq.DEF["RadiofrequencyRasterTime"]
    Δt_gr = seq.DEF["GradientRasterTime"]
    # Get the unique events and its IDs
    rfs = get_unique_events(seq.RF[is_on.(seq.RF)])
    grs = get_unique_events(seq.GR[is_on.(seq.GR)])
    adcs = get_unique_events(seq.ADC[is_on.(seq.ADC)])
    grads, traps = get_unique_gradients(grs)
    rfmags, rfangs, rftimes, cnt = get_rf_shapes(rfs, 1, Δt_rf)
    gramps, grtimes, _ = get_grad_shapes(grads, cnt, Δt_gr)
    # Return the pulseq object
    return (
        blockEvents = get_block_events(seq, rfs, grs, adcs),
        rfLibrary = get_rf_library(rfs, rfmags, rfangs, rftimes, Δt_rf),
        gradLibrary = (
            arb = get_grad_library_arb(grads, gramps, grtimes, Δt_gr),
            trap = get_grad_library_trap(traps),
        ),
        adcLibrary = get_adc_library(adcs),
        shapeLibrary = get_shape_library(rfmags, rfangs, rftimes, gramps, grtimes),
    )
end

"""
    write_seq(seq::Sequence, filename::String)

Writes a .seq file for a given sequence `seq` y the location `filename`
"""
function write_seq(seq::Sequence, filename)
    # Just a warning message for extensions
    #@warn "EXTENSIONS will not be handled"
    # Get the pulseq object
    obj = get_pulseq_object(seq)
    # Write the .seq file
    open(filename, "w") do fid
        @printf(
            fid,
            """
            # Pulseq sequence file
            # Created with KomaMRIFiles.jl %s

            [VERSION]
            major 1
            minor 4
            revision 2
            """,
            string(KomaMRIFiles.__VERSION__)
        )
        if !isempty(seq.DEF)
            @printf(
                fid,
                """

                [DEFINITIONS]
                """
            )
            sorted_keys = sort(collect(keys(seq.DEF)))
            filter!(x -> x != "extension" && x != "additional_text", sorted_keys) # For now remove extensions
            for key in sorted_keys
                val = seq.DEF[key]
                @printf(fid, "%s ", key)
                if isa(val, String)
                    @printf(fid, "%s ", val)
                else
                    if isa(val, Vector{<:Number})
                        for v in val
                            @printf(fid, "%.9g ", v)
                        end
                    else
                        @printf(fid, "%.9g ", val)
                    end
                end
                @printf(fid, "\n")
            end
        end
        if !isempty(obj.blockEvents)
            @printf(
                fid,
                """

                # Format of blocks:
                # NUM DUR RF  GX  GY  GZ  ADC  EXT
                [BLOCKS]
                """
            )
            id_format_str = "%" * string(length(string(length(seq)))) * "d "
            fmt = Printf.Format(id_format_str * "%3d %3d %3d %3d %3d %2d %2d\n")
            for row in obj.blockEvents
                Printf.format(fid, fmt, row...)
            end
        end
        if !isempty(obj.rfLibrary)
            @printf(
                fid,
                """

                # Format of RF events:
                # id amplitude mag_id phase_id time_shape_id delay freq phase
                # ..        Hz   ....     ....          ....    us   Hz   rad
                [RF]
                """
            )
            fmt = Printf.Format("%d %12g %d %d %d %g %g %g\n")
            for row in obj.rfLibrary
                Printf.format(fid, fmt, row...)
            end
        end
        if !isempty(obj.gradLibrary.arb)
            @printf(
                fid,
                """

                # Format of arbitrary gradients:
                #   time_shape_id of 0 means default timing (stepping with grad_raster starting at 1/2 of grad_raster)
                # id amplitude amp_shape_id time_shape_id delay
                # ..      Hz/m       ..         ..          us
                [GRADIENTS]
                """
            )
            for row in obj.gradLibrary.arb
                @printf(fid, "%d %12g %d %d %d\n", row...)
            end
        end
        if !isempty(obj.gradLibrary.trap)
            @printf(
                fid,
                """

                # Format of trapezoid gradients:
                # id amplitude rise flat fall delay
                # ..      Hz/m   us   us   us    us
                [TRAP]
                """
            )
            for row in obj.gradLibrary.trap
                @printf(fid, "%2d %12g %3d %4d %3d %3d\n", row...)
            end
        end
        if !isempty(obj.adcLibrary)
            @printf(
                fid,
                """

                # Format of ADC events:
                # id num dwell delay freq phase
                # ..  ..    ns    us   Hz   rad
                [ADC]
                """
            )
            for row in obj.adcLibrary
                @printf(fid, "%d %d %.0f %.0f %g %g\n", row...)
            end
        end
        if !isempty(obj.shapeLibrary)
            @printf(
                fid,
                """

                # Sequence Shapes
                [SHAPES]

                """
            )
            for row in obj.shapeLibrary
                @printf(fid, "shape_id %d\n", row.id)
                @printf(fid, "num_samples %d\n", row.num)
                [@printf(fid, "%.9g\n", i) for i in row.data]
                @printf(fid, "\n")
            end
        end
    end
    md5hash = bytes2hex(open(md5, filename))
    open(filename, "a") do fid
        @printf(
            fid,
            """

            [SIGNATURE]
            # This is the hash of the Pulseq file, calculated right before the [SIGNATURE] section was added
            # It can be reproduced/verified with md5sum if the file trimmed to the position right above [SIGNATURE]
            # The new line character preceding [SIGNATURE] BELONGS to the signature (and needs to be sripped away for recalculating/verification)
            Type md5
            Hash %s
            """,
            md5hash
        )
    end
end
