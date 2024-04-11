"""
    typeon_obj = get_typeon_obj(seq::Sequence, type::String)

Returns the vector `typeon_obj` of objects (RF, Grad or ADC) where the object is on in a
given sequence `seq` according to the `type` ∈ ["rf", "gr", "adc"] argument.
"""
function get_typeon_obj(seq::Sequence, type::String)
    type == "rf" && return [s.RF[1] for s in seq if is_RF_on(s)]
    type == "gr" && return [g for s in seq for g in s.GR[:, 1] if sum(abs.(g.A)) > 0]
    type == "adc" && return [s.ADC[1] for s in seq if is_ADC_on(s)]
    return []
end

"""
    typeunique_obj_id = get_typeunique_obj_id(typeon_obj::Vector)

Returns the vector `typeunique_obj_id` of objects (RF, Grad or ADC) which are unique given a
input vector `typeon_obj` of objects which are on.
"""
function get_typeunique_obj_id(typeon_obj::Vector)
    typeunique_obj_id = []
    id_cnt = 1
    for obj in typeon_obj
        if all([!(obj ≈ obj_unique) for (obj_unique, _) in typeunique_obj_id])
            push!(typeunique_obj_id, [obj, id_cnt])
            id_cnt += 1
        end
    end
    return typeunique_obj_id
end

"""
    obj = get_obj(s::Sequence, type::String)

Returns the object (RF, Grad or ADC) of a Sequence `s` (ideally a one-block sequence)
according to the input `type` ∈ ["rf", "gx", "gy", "gz", "adc"]).
"""
function get_obj(s::Sequence, type::String)
    type == "rf"  && return s.RF[1]
    type == "gx"  && return s.GR[1,1]
    type == "gy"  && return s.GR[2,1]
    type == "gz"  && return s.GR[3,1]
    type == "adc" && return s.ADC[1]
    return Sequence()
end

"""
    type_blk_obj_id = get_type_blk_obj_id(seq::Sequence, type::String, unique_obj_id::Vector)

Returns the vector `type_blk_obj_id` of vectors [blk, obj, id] for all the blocks `blk` of a
sequence `seq` according to the `type` ∈ ["rf", "gx", "gy", "gz", "adc"] chosen. It is
neccessary to add the input vector `unique_obj_id` which contains the uniques objects (RF,
Grad, or ADC) with its repective ID.
"""
function get_type_blk_obj_id(seq::Sequence, type::String, unique_obj_id::Vector)

    type_blk_obj_id = [[blk, get_obj(s, type), 0] for (blk, s) ∈ enumerate(seq)]
    for boi ∈ type_blk_obj_id
        for (obj_unique, id_unique) ∈ unique_obj_id
            if boi[2] ≈ obj_unique
                boi[3] = id_unique
            end
        end
    end
    return type_blk_obj_id
end

"""
    output = magsign(v)

Returns the extreme value (positive or negative) which has the maximum absolute value.
"""
function magsign(v)
    a, b = extrema(v)
    return -a > b ? a : b
end

"""
    rfunique_abs_id, rfunique_ang_id, rfunique_tim_id, id_shape_cnt = get_rfunique(rfunique_obj_id::Vector, id_shape_cnt::Integer, seq::Sequence)

Returns the unique shapes for the magnitude, angle and time of the "rfunique_obj_id" vector.
Requires an initial integer counter "id_shape_cnt" to asign IDs incrementally.
"""
function get_rfunique(rfunique_obj_id::Vector, id_shape_cnt::Integer, seq::Sequence)
    # Find the unique shapes (magnitude, phase and time shapes) and assign IDs
    rfunique_abs_id, rfunique_ang_id, rfunique_tim_id = [], [], []
    for (obj,_) ∈ rfunique_obj_id
        shape_abs = abs.(obj.A) / maximum(abs.(obj.A))
        if all([!(length(shape_abs) == length(shape_abs_unique) && shape_abs ≈ shape_abs_unique) for (shape_abs_unique,_) ∈ rfunique_abs_id])
            push!(rfunique_abs_id, [shape_abs, id_shape_cnt])
            id_shape_cnt += 1
        end
        shape_ang = mod.(angle.(obj.A), 2π)/2π
        if all([!(length(shape_ang) == length(shape_ang_unique) && all(shape_ang - shape_ang_unique .≈ shape_ang[1] - shape_ang_unique[1])) for (shape_ang_unique,_) ∈ rfunique_ang_id])
            push!(rfunique_ang_id, [shape_ang, id_shape_cnt])
            id_shape_cnt += 1
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / seq.DEF["RadiofrequencyRasterTime"]
            if all([!(length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique) for (shape_tim_unique,_) ∈ rfunique_tim_id])
                push!(rfunique_tim_id, [shape_tim, id_shape_cnt])
                id_shape_cnt += 1
            end
        end
    end
    return rfunique_abs_id, rfunique_ang_id, rfunique_tim_id, id_shape_cnt
end


"""
    gradunique_amp_id, gradunique_tim_id, id_shape_cnt = get_gradunique(gradunique_obj_id::Vector, id_shape_cnt::Integer, seq::Sequence)

Returns the unique shapes for the amplitude and time of the "gradunique_obj_id" vector.
Requires an initial integer counter "id_shape_cnt" to asign IDs incrementally.
"""
function get_gradunique(gradunique_obj_id::Vector, id_shape_cnt::Integer, seq::Sequence)
    # Find shapes for magnitude and time gradients
    gradunique_amp_id, gradunique_tim_id = [], []
    for (obj,_) ∈ gradunique_obj_id
        shape_amp = obj.A / magsign(obj.A)
        if all([!(length(shape_amp) == length(shape_amp_unique) && shape_amp ≈ shape_amp_unique) for (shape_amp_unique,_) ∈ gradunique_amp_id])
            push!(gradunique_amp_id, [shape_amp, id_shape_cnt])
            id_shape_cnt = id_shape_cnt + 1
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T])/seq.DEF["GradientRasterTime"]
            if all([!(length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique) for (shape_tim_unique,_) ∈ gradunique_tim_id])
                push!(gradunique_tim_id, [shape_tim, id_shape_cnt])
                id_shape_cnt = id_shape_cnt + 1
            end
        end
    end
    return gradunique_amp_id, gradunique_tim_id, id_shape_cnt
end

"""
    compressed_data = compress(data::Vector)

Returns the compressed data vector `compressed_data` for the intupt vector `data` using the
algorithm for encoding the derivative in a run-length compressed format.
"""
function compress(data)
    return compress_shape(data)[2]
end

"""
    write_seq(seq::Sequence, filename::String)

Writes a .seq file for a given sequence `seq` y the location `filename`
"""
function write_seq(seq::Sequence, filename)

    # Get the unique objects (RF, Grad y ADC) and its IDs
    rfunique_obj_id  = get_typeunique_obj_id(get_typeon_obj(seq, "rf"))
    grunique_obj_id  = get_typeunique_obj_id(get_typeon_obj(seq, "gr"))
    adcunique_obj_id = get_typeunique_obj_id(get_typeon_obj(seq, "adc"))
    gradunique_obj_id = [[obj, id] for (obj, id) ∈ grunique_obj_id if length(obj.A) != 1]
    trapunique_obj_id = [[obj, id] for (obj, id) ∈ grunique_obj_id if length(obj.A) == 1]
    rfunique_abs_id, rfunique_ang_id, rfunique_tim_id, id_shape_cnt = get_rfunique(rfunique_obj_id, 1, seq)
    gradunique_amp_id, gradunique_tim_id, _ = get_gradunique(gradunique_obj_id, id_shape_cnt, seq)

    # Define the table to be written for the [BLOCKS] section
    @warn "EXTENSIONS will not be handled"
    rf_id  = [id for (_, _, id) ∈ get_type_blk_obj_id(seq, "rf",  rfunique_obj_id)]
    gx_id  = [id for (_, _, id) ∈ get_type_blk_obj_id(seq, "gx",  grunique_obj_id)]
    gy_id  = [id for (_, _, id) ∈ get_type_blk_obj_id(seq, "gy",  grunique_obj_id)]
    gz_id  = [id for (_, _, id) ∈ get_type_blk_obj_id(seq, "gz",  grunique_obj_id)]
    adc_id = [id for (_, _, id) ∈ get_type_blk_obj_id(seq, "adc", adcunique_obj_id)]
    blk_obj_dur_rf_gx_gy_gz_adc_ext = [[b, s, 0, rf_id[b], gx_id[b], gy_id[b], gz_id[b], adc_id[b], 0] for (b, s) ∈ enumerate(seq)]
    for bodrxyzae ∈ blk_obj_dur_rf_gx_gy_gz_adc_ext
        blk = bodrxyzae[1]
        bd = seq.DUR[blk] / seq.DEF["BlockDurationRaster"]
        bdr = round(bd)
        if abs(bdr-bd) > 1e-6
            @warn "Block $blk duration rounded"
        end
        bodrxyzae[3] = bdr
    end

    # Define the table to be written for the [RF] section
    rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha = [[idx, obj, 0, 0, 0, 0, 0, 0, 0] for (obj, idx) ∈ rfunique_obj_id]
    for ioamptdfh ∈ rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha
        obj = ioamptdfh[2]
        ioamptdfh[3] = γ * maximum(abs.(obj.A))
        ioamptdfh[8] = obj.Δf
        shape_abs = abs.(obj.A) / maximum(abs.(obj.A))
        for (shape_abs_unique, id_abs) ∈ rfunique_abs_id
            if length(shape_abs) == length(shape_abs_unique) && shape_abs ≈ shape_abs_unique
                ioamptdfh[4] = id_abs
            end
        end
        shape_ang = mod.(angle.(obj.A), 2π)/2π
        for (shape_ang_unique, id_ang) ∈ rfunique_ang_id
            if length(shape_ang) == length(shape_ang_unique) && all(shape_ang - shape_ang_unique .≈ shape_ang[1] - shape_ang_unique[1])
                ioamptdfh[5] = id_ang
                ioamptdfh[9] = angle(sum(exp.(1im*2π*shape_ang) .* exp.(-1im*2π*shape_ang_unique))/length(shape_ang))
            end
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T])/seq.DEF["RadiofrequencyRasterTime"]
            for (shape_tim_unique, id_tim) ∈ rfunique_tim_id
                if length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique
                    ioamptdfh[6] = id_tim
                end
            end
        end
        delay_compensation_rf_koma = (ioamptdfh[6] == 0) * seq.DEF["RadiofrequencyRasterTime"] / 2
        ioamptdfh[7] = round((obj.delay - delay_compensation_rf_koma) / seq.DEF["RadiofrequencyRasterTime"]) * seq.DEF["RadiofrequencyRasterTime"] * 1e6
    end

    # Define the table to be written for the [GRADIENTS] section
    grad_idx_obj_amp_iamp_itim_delay = [
        [idx, obj, 0, 0, 0, 0] for (obj, idx) in gradunique_obj_id
    ]
    for ioamtd in grad_idx_obj_amp_iamp_itim_delay
        obj = ioamtd[2]
        ioamtd[3] = γ * magsign(obj.A)
        ioamtd[6] = round(1e6 * obj.delay)
        shape_amp = obj.A / magsign(obj.A)
        for (shape_amp_unique, id_amp) in gradunique_amp_id
            if length(shape_amp) == length(shape_amp_unique) && shape_amp ≈ shape_amp_unique
                ioamtd[4] = id_amp
            end
        end
        if isa(obj.T, Vector{<:Number})
            shape_tim = cumsum([0; obj.T]) / seq.DEF["GradientRasterTime"]
            for (shape_tim_unique, id_tim) ∈ gradunique_tim_id
                if length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique
                    ioamtd[5] = id_tim
                end
            end
        end
    end

    # Define the table to be written for the [TRAP] section
    trap_idx_obj_amp_rise_flat_fall_delay = [[idx, obj, 0, 0, 0, 0, 0] for (obj, idx) ∈ trapunique_obj_id]
    for ioarfad ∈ trap_idx_obj_amp_rise_flat_fall_delay
        obj = ioarfad[2]
        ioarfad[3] = γ * obj.A
        ioarfad[4] = 1e6 * obj.rise
        ioarfad[5] = 1e6 * obj.T
        ioarfad[6] = 1e6 * obj.fall
        ioarfad[7] = 1e6 * obj.delay
    end

    # Define the table to be written for the [ADC] section
    adc_idx_obj_num_dwell_delay_freq_phase = [[idx, obj, 0, 0, 0, 0, 0] for (obj, idx) ∈ adcunique_obj_id]
    for ionwdfp ∈ adc_idx_obj_num_dwell_delay_freq_phase
        obj = ionwdfp[2]
        ionwdfp[3] = obj.N
        ionwdfp[4] = obj.T * 1e9 / (obj.N - 1)
        ionwdfp[5] = (obj.delay - 0.5*obj.T/(obj.N - 1)) * 1e6
        ionwdfp[6] = obj.Δf
        ionwdfp[7] = obj.ϕ
    end

    # Define the table to be written for the [SHAPES] section
    shapefull_data_id = [shapeunique_data_id_i for shapeunique_data_id ∈ [rfunique_abs_id, rfunique_ang_id, rfunique_tim_id, gradunique_amp_id, gradunique_tim_id] for shapeunique_data_id_i ∈ shapeunique_data_id]
    shape_data_id_num = [(length(compress(data)) == length(data) ? data : compress(data), id, length(data)) for (data, id) ∈ shapefull_data_id]

    # Write the .seq file
    open(filename, "w") do fid

        @printf(fid, "# Pulseq sequence file\n")
        @printf(fid, "# Created by KomaMRI.jl \n\n") #TODO: add Koma version

        @printf(fid, "[VERSION]\n")
        @printf(fid, "major 1\n")
        @printf(fid, "minor 4\n")
        @printf(fid, "revision 1\n")
        @printf(fid, "\n")

        if !isempty(seq.DEF)
            @printf(fid, "[DEFINITIONS]\n")
            sorted_keys = sort(collect(keys(seq.DEF)))
            for key ∈ sorted_keys
                val = seq.DEF[key]
                if key ∈ ["AdcRasterTime", "BlockDurationRaster", "GradientRasterTime", "RadiofrequencyRasterTime", "TotalDuration", "FOV", "MaxAdcSegmentLength"]
                    @printf(fid, "%s ", key)
                    if isa(val, String)
                        @printf(fid, "%s ", val)
                    else
                        if isa(val, Vector{<:Number})
                            for v ∈ val
                                @printf(fid, "%.9g ", v)
                            end
                        else
                            @printf(fid, "%.9g ", val)
                        end
                    end
                    @printf(fid, "\n")
                end
            end
            @printf(fid, "\n")
        end

        if !isempty(blk_obj_dur_rf_gx_gy_gz_adc_ext)
            @printf(fid, "# Format of blocks:\n")
            @printf(fid, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
            @printf(fid, "[BLOCKS]\n")
            nBlocks = length(seq)
            idFormatWidth = length(string(nBlocks))
            idFormatStr = "%" * string(idFormatWidth) * "d "
            for (blk, _, dur, rf, gx, gy, gz, adc, ext) ∈ blk_obj_dur_rf_gx_gy_gz_adc_ext
                Printf.format(fid, Printf.Format(idFormatStr * "%3d %3d %3d %3d %3d %2d %2d\n"), blk, dur, rf, gx, gy, gz, adc, ext)
            end
            @printf(fid, "\n")
        end

        if !isempty(rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha)
            @printf(fid, "# Format of RF events:\n")
            @printf(fid, "# id amplitude mag_id phase_id time_shape_id delay freq phase\n")
            @printf(fid, "# ..        Hz   ....     ....          ....    us   Hz   rad\n")
            @printf(fid, "[RF]\n")
            for (id, _, amp, magid, phaid, timeid, delay, freq, phase) ∈ rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha
                @printf(fid, "%d %12g %d %d %d %g %g %g\n", id, amp, magid, phaid, timeid, delay, freq, phase)
            end
            @printf(fid, "\n")
        end

        if !isempty(grad_idx_obj_amp_iamp_itim_delay)
            @printf(fid, "# Format of arbitrary gradients:\n")
            @printf(fid, "#   time_shape_id of 0 means default timing (stepping with grad_raster starting at 1/2 of grad_raster)\n")
            @printf(fid, "# id amplitude amp_shape_id time_shape_id delay\n") # do we need delay ???
            @printf(fid, "# ..      Hz/m       ..         ..          us\n")
            @printf(fid, "[GRADIENTS]\n")
            for (id, _, amp, ampid, timeid, delay) ∈ grad_idx_obj_amp_iamp_itim_delay
                @printf(fid, "%d %12g %d %d %d\n", id, amp, ampid, timeid, delay)
            end
            @printf(fid, "\n")
        end

        if !isempty(trap_idx_obj_amp_rise_flat_fall_delay)
            @printf(fid, "# Format of trapezoid gradients:\n")
            @printf(fid, "# id amplitude rise flat fall delay\n")
            @printf(fid, "# ..      Hz/m   us   us   us    us\n")
            @printf(fid, "[TRAP]\n")
            for (id, _, amp, rise, flat, fall, delay) ∈ trap_idx_obj_amp_rise_flat_fall_delay
                @printf(fid, "%2d %12g %3d %4d %3d %3d\n", id, amp, rise, flat, fall, delay)
            end
            @printf(fid, "\n")
        end

        if !isempty(adc_idx_obj_num_dwell_delay_freq_phase)
            @printf(fid, "# Format of ADC events:\n")
            @printf(fid, "# id num dwell delay freq phase\n")
            @printf(fid, "# ..  ..    ns    us   Hz   rad\n")
            @printf(fid, "[ADC]\n")
            for (id, _, num, dwell, delay, freq, phase) ∈ adc_idx_obj_num_dwell_delay_freq_phase
                @printf(fid, "%d %d %.0f %.0f %g %g\n", id, num, dwell, delay, freq, phase)
            end
            @printf(fid, "\n")
        end

        if !isempty(shape_data_id_num)
            @printf(fid, "# Sequence Shapes\n")
            @printf(fid, "[SHAPES]\n\n")
            for (data, id, num) ∈ shape_data_id_num
                @printf(fid, "shape_id %d\n", id)
                @printf(fid, "num_samples %d\n", num)
                [@printf(fid, "%.9g\n", datai) for datai ∈ data]
                @printf(fid, "\n")
            end
        end

    end

    md5hash = bytes2hex(open(md5, filename))
    open(filename, "a") do fid
        @printf(fid, "\n[SIGNATURE]\n") # the preceding new line BELONGS to the signature (and needs to be sripped away to recalculate the signature)
        @printf(fid, "# This is the hash of the Pulseq file, calculated right before the [SIGNATURE] section was added\n")
        @printf(fid, "# It can be reproduced/verified with md5sum if the file trimmed to the position right above [SIGNATURE]\n")
        @printf(fid, "# The new line character preceding [SIGNATURE] BELONGS to the signature (and needs to be sripped away for recalculating/verification)\n")
        @printf(fid, "Type md5\n")
        @printf(fid, "Hash %s\n", md5hash)
    end

end
