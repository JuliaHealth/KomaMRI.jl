using KomaMRI, Printf, MD5

function magsign(v)
    a, b = extrema(v)
    return -a > b ? a : b
end

function compress(data)
    derivative = [[data[1]]; data[2:end] - data[1:end-1]]
    encoded = Tuple{Float64, Int64}[]
    current_value = derivative[1]
    count = 1
    for i in 2:length(derivative)
        if derivative[i] == current_value
            count += 1
        else
            push!(encoded, (current_value, count))
            current_value = derivative[i]
            count = 1
        end
    end
    push!(encoded, (current_value, count))

    compression = []
    for (v,r) in encoded
        if r == 1
            push!(compression, v)
        elseif r == 2
            push!(compression, v)
            push!(compression, v)
        elseif r > 2
            push!(compression, v)
            push!(compression, v)
            push!(compression, Int(r-2))
        end
    end

    return compression
end


Base.isapprox(rf1::RF, rf2::RF) = begin
    return all(length(getfield(rf1,k)) == length(getfield(rf2,k)) for k ∈ fieldnames(RF)) &&
        all(size(getfield(rf1,k)) == size(getfield(rf2,k)) for k ∈ fieldnames(RF)) &&
        all(≈(getfield(rf1,k), getfield(rf2,k), atol=1e-9) for k ∈ fieldnames(RF))
end

Base.isapprox(gr1::Grad, gr2::Grad) = begin
    return all(length(getfield(gr1,k)) ≈ length(getfield(gr2,k)) for k ∈ fieldnames(Grad)) &&
        #all(size(getfield(gr1,k)) ≈ size(getfield(gr2,k)) for k ∈ fieldnames(Grad)) &&
        all(getfield(gr1,k) ≈ getfield(gr2,k) for k ∈ fieldnames(Grad))
end

Base.isapprox(adc1::ADC, adc2::ADC) = begin
    return all(length(getfield(adc1,k)) ≈ length(getfield(adc2,k)) for k ∈ fieldnames(ADC)) &&
        #all(size(getfield(adc1,k)) ≈ size(getfield(adc2,k)) for k ∈ fieldnames(ADC)) &&
        all(getfield(adc1,k) ≈ getfield(adc2,k) for k ∈ fieldnames(ADC))
end

Base.isapprox(s1::Sequence, s2::Sequence) = begin
    if length(s1) != length(s2)
        return false
    end
    return all([s1.ADC[i] ≈ s2.ADC[i] && s1.RF[i] ≈ s2.RF[i] && s1.GR[i] ≈ s2.GR[i] && s1.DUR[i] ≈ s2.DUR[i] for i in 1:length(s1)])
end

id_cnt = 1
id_shape_cnt = 1

function write_seq(seq::Sequence, filename)

# Find the objects where the rf is on
rfon_obj = [s.RF[1] for s ∈ seq if is_RF_on(s)]
# Find the unique rfs and assign IDs
global id_cnt = 1
rfunique_obj_id = []
for obj ∈ rfon_obj
    global id_cnt
    if all([!(obj ≈ obj_unique) for (obj_unique,_) ∈ rfunique_obj_id])
        push!(rfunique_obj_id, [obj, id_cnt])
        id_cnt = id_cnt + 1
    end
end
# Fill the IDs for all the rfs (for [BLOCKS] section)
rf_blk_obj_id = [[b, s.RF[1], 0] for (b,s) ∈ enumerate(seq)]
for boi ∈ rf_blk_obj_id
    for (obj_unique, id_unique) ∈ rfunique_obj_id
        if boi[2] ≈ obj_unique
            boi[3] = id_unique
        end
    end
end

# Find the shape (magnitude, phase and time shapes) of the rfs when they are on and possibly with unique shapes
rfon_abs = [abs.(obj.A)/maximum(abs.(obj.A)) for (obj,_) ∈ rfunique_obj_id]
rfon_ang = [-angle.(obj.A)/2π for (obj,_) ∈ rfunique_obj_id]
rfon_tim = [cumsum([0; obj.T])/seq.DEF["RadiofrequencyRasterTime"] for (obj,_) ∈ rfunique_obj_id if isa(obj.T, Vector{<:Number})]
# Find the unique shapes (magnitude, phase and time shapes) and assign IDs
global id_shape_cnt = 1
rfunique_abs_id, rfunique_ang_id, rfunique_tim_id = [], [], []
for (obj,_) ∈ rfunique_obj_id
    global id_shape_cnt
    shape_abs = abs.(obj.A)/maximum(abs.(obj.A))
    if all([!(length(shape_abs) == length(shape_abs_unique) && shape_abs ≈ shape_abs_unique) for (shape_abs_unique,_) ∈ rfunique_abs_id])
        push!(rfunique_abs_id, [shape_abs, id_shape_cnt])
        id_shape_cnt = id_shape_cnt + 1
    end
    shape_ang = -angle.(obj.A)/2π
    if all([!(length(shape_ang) == length(shape_ang_unique) && all(shape_ang - shape_ang_unique .≈ shape_ang[1] - shape_ang_unique[1])) for (shape_ang_unique,_) ∈ rfunique_ang_id])
        push!(rfunique_ang_id, [shape_ang, id_shape_cnt])
        id_shape_cnt = id_shape_cnt + 1
    end
    if isa(obj.T, Vector{<:Number})
        shape_tim = cumsum([0; obj.T])/seq.DEF["RadiofrequencyRasterTime"]
        if all([!(length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique) for (shape_tim_unique,_) ∈ rfunique_tim_id])
            push!(rfunique_tim_id, [shape_tim, id_shape_cnt])
            id_shape_cnt = id_shape_cnt + 1
        end
    end
end
# Fill the IDs for all the rf shapes and some additional values (for [RF] section)
rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha = [[idx, obj, 0, 0, 0, 0, 0, 0, 0] for (obj, idx) ∈ rfunique_obj_id]
for ioamptdfh ∈ rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha
    obj = ioamptdfh[2]
    ioamptdfh[3] = γ * maximum(abs.(obj.A))
    ioamptdfh[7] = round(obj.delay / seq.DEF["RadiofrequencyRasterTime"]) * seq.DEF["RadiofrequencyRasterTime"] * 1e6
    ioamptdfh[8] = obj.Δf
    shape_abs = abs.(obj.A)/maximum(abs.(obj.A))
    for (shape_abs_unique, id_abs) ∈ rfunique_abs_id
        if length(shape_abs) == length(shape_abs_unique) && shape_abs ≈ shape_abs_unique
            ioamptdfh[4] = id_abs
        end
    end
    shape_ang = -angle.(obj.A)/2π
    for (shape_ang_unique, id_ang) ∈ rfunique_ang_id
        if length(shape_ang) == length(shape_ang_unique) && all(shape_ang - shape_ang_unique .≈ shape_ang[1] - shape_ang_unique[1])
            ioamptdfh[5] = id_ang
            a = 2π*(shape_ang[1] + shape_ang_unique[1])
            if a < 0
                a += 2π
            end
            ioamptdfh[9] = a
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
end

### For GR
# Find the objects where the gr is on
gron_obj = [g for s in seq for g in s.GR[:,1] if sum(abs.(g.A)) > 0]
# Find the unique grs and assign IDs
id_cnt = 1
grunique_obj_id = []
for obj ∈ gron_obj
    global id_cnt
    if all([!(obj ≈ obj_unique) for (obj_unique,_) ∈ grunique_obj_id])
        push!(grunique_obj_id, [obj, id_cnt])
        id_cnt = id_cnt + 1
    end
end
### For GX:
gxon_blk_obj = [[blk, s.GR[1,1]] for (blk, s) in enumerate(seq) if KomaMRICore.is_Gx_on(s)]
gxon_blk_obj_id = [[blk, obj, id] for (blk, obj) in gxon_blk_obj for (obj_unique, id) in grunique_obj_id if obj ≈ obj_unique]
gx_blk_obj_id = [[blk, s.GR[1,1], 0] for (blk,s) in enumerate(seq)]
[boi[3] = id for boi in gx_blk_obj_id for (blk, _, id) in gxon_blk_obj_id if boi[1] == blk]
### For GY:
gyon_blk_obj = [[blk, s.GR[2,1]] for (blk, s) in enumerate(seq) if KomaMRICore.is_Gy_on(s)]
gyon_blk_obj_id = [[blk, obj, id] for (blk, obj) in gyon_blk_obj for (obj_unique, id) in grunique_obj_id if obj ≈ obj_unique]
gy_blk_obj_id = [[blk, s.GR[2,1], 0] for (blk,s) in enumerate(seq)]
[boi[3] = id for boi in gy_blk_obj_id for (blk, _, id) in gyon_blk_obj_id if boi[1] == blk]
### For GZ:
gzon_blk_obj = [[blk, s.GR[3,1]] for (blk, s) in enumerate(seq) if KomaMRICore.is_Gz_on(s)]
gzon_blk_obj_id = [[blk, obj, id] for (blk, obj) in gzon_blk_obj for (obj_unique, id) in grunique_obj_id if obj ≈ obj_unique]
gz_blk_obj_id = [[blk, s.GR[3,1], 0] for (blk,s) in enumerate(seq)]
[boi[3] = id for boi in gz_blk_obj_id for (blk, _, id) in gzon_blk_obj_id if boi[1] == blk]


### For GRAD
gradunique_obj_id = [[obj, id] for (obj, id) in grunique_obj_id if length(obj.A) != 1]
gradon_amp = [obj.A / magsign(obj.A) for (obj,_) ∈ gradunique_obj_id]
gradon_tim = [cumsum([0; obj.T])/seq.DEF["GradientRasterTime"] for (obj,_) ∈ gradunique_obj_id if isa(obj.T, Vector{<:Number})]
#
gradunique_amp_id, gradunique_tim_id = [], []
for (obj,_) ∈ gradunique_obj_id
    global id_shape_cnt
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
#
grad_idx_obj_amp_iamp_itim_delay = [[idx, obj, 0, 0, 0, 0] for (obj, idx) ∈ gradunique_obj_id]
for ioamtd ∈ grad_idx_obj_amp_iamp_itim_delay
    obj = ioamtd[2]
    ioamtd[3] = γ * magsign(obj.A)
    ioamtd[6] = round(1e6 * obj.delay)
    shape_amp = obj.A / magsign(obj.A)
    for (shape_amp_unique, id_amp) ∈ gradunique_amp_id
        if length(shape_amp) == length(shape_amp_unique) && shape_amp ≈ shape_amp_unique
            ioamtd[4] = id_amp
        end
    end
    if isa(obj.T, Vector{<:Number})
        shape_tim = cumsum([0; obj.T])/seq.DEF["GradientRasterTime"]
        for (shape_tim_unique, id_tim) ∈ gradunique_tim_id
            if length(shape_tim) == length(shape_tim_unique) && shape_tim ≈ shape_tim_unique
                ioamtd[5] = id_tim
            end
        end
    end
end

### For TRAP
trap_obj_id = [[obj, id] for (obj, id) in grunique_obj_id if length(obj.A) == 1]
trap_idx_obj_amp_rise_flat_fall_delay = [[idx, obj, 0, 0, 0, 0, 0] for (obj, idx) ∈ trap_obj_id]
for ioarfad ∈ trap_idx_obj_amp_rise_flat_fall_delay
    obj = ioarfad[2]
    ioarfad[3] = γ * obj.A
    ioarfad[4] = 1e6 * obj.rise
    ioarfad[5] = 1e6 * obj.T
    ioarfad[6] = 1e6 * obj.fall
    ioarfad[7] = 1e6 * obj.delay
end

### For ADC
adcon_blk_obj = [(blk, s.ADC[1]) for (blk,s) in enumerate(seq) if is_ADC_on(s)]
id_cnt = 1
adcunique_blk_obj_id = []
for (blk, obj) ∈ adcon_blk_obj
    global id_cnt
    if all([!(obj ≈ obj_unique) for (_, obj_unique, _) ∈ adcunique_blk_obj_id])
        push!(adcunique_blk_obj_id, [blk, obj, id_cnt])
        id_cnt = id_cnt + 1
    end
end
adcon_blk_obj_id = [[blk, obj, id] for (blk, obj) in adcon_blk_obj for (_,obj_unique, id) in adcunique_blk_obj_id if obj ≈ obj_unique]
adc_blk_obj_id = [[blk, s.ADC[1], 0] for (blk,s) in enumerate(seq)]
[boi[3] = id for boi in adc_blk_obj_id for (blk, _, id) in adcon_blk_obj_id if boi[1] == blk]
#
adc_idx_obj_num_dwell_delay_freq_phase = [[idx, obj, 0, 0, 0, 0, 0] for (_, obj, idx) ∈ adcunique_blk_obj_id]
for ionwdfp ∈ adc_idx_obj_num_dwell_delay_freq_phase
    obj = ionwdfp[2]
    ionwdfp[3] = obj.N
    ionwdfp[4] = obj.T * 1e9 / (obj.N - 1)
    ionwdfp[5] = (obj.delay - 0.5*obj.T/(obj.N - 1)) * 1e6
    ionwdfp[6] = obj.Δf
    ionwdfp[7] = obj.ϕ
end

### For SHAPES
# Fix last sample of RF
[push!(a[1], a[1][end]) for a ∈ rfunique_abs_id]
[push!(a[1], a[1][end]) for a ∈ rfunique_ang_id]
shapefull_data_id = []
[push!(shapefull_data_id, a) for a ∈ rfunique_abs_id]
[push!(shapefull_data_id, a) for a ∈ rfunique_ang_id]
[push!(shapefull_data_id, a) for a ∈ rfunique_tim_id]
[push!(shapefull_data_id, a) for a ∈ gradunique_amp_id]
[push!(shapefull_data_id, a) for a ∈ gradunique_tim_id]
shape_data_id_num = [(length(compress(data)) == length(data) ? data : compress(data), id, length(data)) for (data, id) in shapefull_data_id]

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
        for key in sorted_keys
            val = seq.DEF[key]
            if key in ["AdcRasterTime"; "BlockDurationRaster"; "GradientRasterTime"; "RadiofrequencyRasterTime"; "TotalDuration"; "FOV"; "MaxAdcSegmentLength"]
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
                @printf(fid, "\n");
            end
        end
        @printf(fid, "\n")
    end

    @printf(fid, "# Format of blocks:\n")
    @printf(fid, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
    @printf(fid, "[BLOCKS]\n")
    iext = 0
    @warn "EXTENSIONS will not be handled"
    nBlocks = length(seq)
    idFormatWidth = length(string(nBlocks))
    idFormatStr = "%" * string(idFormatWidth) * "d "
    for (i, s) in enumerate(seq)
        bd = seq.DUR[i] / seq.DEF["BlockDurationRaster"]
        bdr = round(bd)
        if abs(bdr-bd) > 1e-6
            @warn "Block $i duration rounded"
        end
        irf = rf_blk_obj_id[i][3]
        igx = gx_blk_obj_id[i][3]
        igy = gy_blk_obj_id[i][3]
        igz = gz_blk_obj_id[i][3]
        iadc = adc_blk_obj_id[i][3]
        Printf.format(fid, Printf.Format(idFormatStr * "%3d %3d %3d %3d %3d %2d %2d\n"), i, bdr, irf, igx, igy, igz, iadc, iext)
    end
    @printf(fid, "\n")

    if !isempty(rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha)
        @printf(fid, "# Format of RF events:\n")
        @printf(fid, "# id amplitude mag_id phase_id time_shape_id delay freq phase\n")
        @printf(fid, "# ..        Hz   ....     ....          ....    us   Hz   rad\n")
        @printf(fid, "[RF]\n")
        for (id, _, amp, magid, phaid, timeid, delay, freq, phase) in rf_idx_obj_amp_imag_ipha_itim_delay_freq_pha
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
        for (id, _, amp, ampid, timeid, delay) in grad_idx_obj_amp_iamp_itim_delay
            @printf(fid, "%d %12g %d %d %d\n", id, amp, ampid, timeid, delay)
        end
        @printf(fid, "\n")
    end

    if !isempty(trap_idx_obj_amp_rise_flat_fall_delay)
        @printf(fid, "# Format of trapezoid gradients:\n")
        @printf(fid, "# id amplitude rise flat fall delay\n")
        @printf(fid, "# ..      Hz/m   us   us   us    us\n")
        @printf(fid, "[TRAP]\n")
        for (id, _, amp, rise, flat, fall, delay) in trap_idx_obj_amp_rise_flat_fall_delay
            @printf(fid, "%2d %12g %3d %4d %3d %3d\n", id, amp, rise, flat, fall, delay)
        end
        @printf(fid, "\n")
    end

    if !isempty(adc_idx_obj_num_dwell_delay_freq_phase)
        @printf(fid, "# Format of ADC events:\n")
        @printf(fid, "# id num dwell delay freq phase\n")
        @printf(fid, "# ..  ..    ns    us   Hz   rad\n")
        @printf(fid, "[ADC]\n")
        for (id, _, num, dwell, delay, freq, phase) in adc_idx_obj_num_dwell_delay_freq_phase
            @printf(fid, "%d %d %.0f %.0f %g %g\n", id, num, dwell, delay, freq, phase)
        end
        @printf(fid, "\n")
    end

    if !isempty(shape_data_id_num)
        @printf(fid, "# Sequence Shapes\n")
        @printf(fid, "[SHAPES]\n\n")
        for (data, id, num) in shape_data_id_num
            @printf(fid, "shape_id %d\n", id)
            @printf(fid, "num_samples %d\n", num)
            [@printf(fid, "%.9g\n", datai) for datai in data]
            @printf(fid, "\n")
        end
    end

end

md5hash = bytes2hex(open(md5, filename))
open(filename, "a") do fid
    @printf(fid, "\n[SIGNATURE]\n"); # the preceding new line BELONGS to the signature (and needs to be sripped away to recalculate the signature)
    @printf(fid, "# This is the hash of the Pulseq file, calculated right before the [SIGNATURE] section was added\n");
    @printf(fid, "# It can be reproduced/verified with md5sum if the file trimmed to the position right above [SIGNATURE]\n");
    @printf(fid, "# The new line character preceding [SIGNATURE] BELONGS to the signature (and needs to be sripped away for recalculating/verification)\n");
    @printf(fid, "Type md5\n");
    @printf(fid, "Hash %s\n", md5hash);
end

end

seq_file = "epi0.seq"
#seq_file = "epi_se0.seq"
#seq_file = "spiral0.seq"
#seq_file = "gre_rad0.seq"
#seq_file = "cine_gre0.seq"

seq = read_seq(seq_file)
filename = splitext(seq_file)[1][1:end-1] * "1.seq"
write_seq(seq, filename)
seq1 = read_seq(filename)
println("seq ≈ seq1")
println(seq ≈ seq1)
