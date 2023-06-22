using KomaMRI, Printf, MD5

#seq_file = "epi0.seq"
#seq_file = "spiral0.seq"
seq_file = "gre_rad0.seq"
seq = read_seq(seq_file)

plot_seq(seq)
plot_kspace(seq)

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
    return all(length(getfield(rf1,k)) ≈ length(getfield(rf2,k)) for k ∈ fieldnames(RF)) &&
        all(getfield(rf1,k) ≈ getfield(rf2,k) for k ∈ fieldnames(RF))
end

Base.isapprox(gr1::Grad, gr2::Grad) = begin
    return all(length(getfield(gr1,k)) ≈ length(getfield(gr2,k)) for k ∈ fieldnames(Grad)) &&
        all(getfield(gr1,k) ≈ getfield(gr2,k) for k ∈ fieldnames(Grad))
end

Base.isapprox(adc1::ADC, adc2::ADC) = begin
    return all(length(getfield(adc1,k)) ≈ length(getfield(adc2,k)) for k ∈ fieldnames(ADC)) &&
        all(getfield(adc1,k) ≈ getfield(adc2,k) for k ∈ fieldnames(ADC))
end


filename = splitext(seq_file)[1][1:end-1] * "1.seq"
open(filename, "w") do fid

    @printf(fid, "# Pulseq sequence file\n")
    @printf(fid, "# Created by KomaMRI\n\n")

    @printf(fid, "[VERSION]\n")
    @printf(fid, "major 1\n")
    @printf(fid, "minor 4\n")
    @printf(fid, "revision 1\n")
    @printf(fid, "\n")

    if !isempty(seq.DEF)
        @printf(fid, "[DEFINITIONS]\n")
        for (key, val) in seq.DEF
            if key in ["AdcRasterTime"; "BlockDurationRaster"; "GradientRasterTime"; "RadiofrequencyRasterTime"; "TotalDuration"]
                @printf(fid, "%s ", key)
                if isa(val, String)
                    @printf(fid, "%s ", val)
                else
                    @printf(fid, "%.9g ", val)
                end
                @printf(fid, "\n");
            end
        end
        @printf(fid, "\n")
    end

    # For RF:
    borfon = [(b, s.RF[1]) for (b,s) in enumerate(seq) if is_RF_on(s)]
    boirfun = [(borfon[1][1], borfon[1][2], 1)]
    [push!(boirfun, (b, o, boirfun[end][3]+1)) for (b,o) in borfon if all([!(o ≈ ou) for (bu,ou,iu) in boirfun])]
    boirfon = [(b,o,iu) for (b,o) in borfon for (bu,ou,iu) in boirfun if o ≈ ou]
    boirf = [[b, s.RF[1], 0] for (b,s) in enumerate(seq)]
    [boi[3] = io for boi in boirf for (bo,oo,io) in boirfon if boi[1] == bo]
    #
    oiawvrfun = [(o, i, maximum(abs.(o.A))*γ, o.A / maximum(abs.(o.A)), cumsum([0; o.T])/seq.DEF["RadiofrequencyRasterTime"]) for (_,o,i) in boirfun]
    wishapeun = [(oiawvrfun[1][4], 1)]
    #[push!(wishapeun, (w, wishapeun[end][2]+1)) for (o,i,a,w,v) in oiawvrfun if all([!(length(w) == length(ws) && w ≈ ws) for (ws,is) in wishapeun])]
    #[push!(wishapeun, (w, wishapeun[end][2]+1)) for (o,i,a,w,v) in oiawvrfun if all([!(length(w) == length(ws) && abs.(w) ≈ abs.(ws) && (all(angle.(w) - angle.(ws) .≈ angle(w[1]) - angle.(ws[1])))) for (ws,is) in wishapeun])]
    [push!(wishapeun, (w, wishapeun[end][2]+1)) for (o,i,a,w,v) in oiawvrfun if all([!(length(w) == length(ws) && abs.(w) ≈ abs.(ws) && (all(angle.(w) - angle.(ws) .≈ angle(w[1]) - angle.(ws[1])))) for (ws,is) in wishapeun])]
    # Assign IDs for RF shapes (this assumes always 2 IDs for magnitude and phase for every rf shape):
    wishape = [(abs.(wishapeun[1][1]), 1); (-angle.(wishapeun[1][1])/2π, 2)]
    [push!(wishape, (abs.(w), wishape[end][2]+1), (-angle.(w)/2π, wishape[end][2]+2)) for (k,(w,i)) in enumerate(wishapeun) if k != 1]
    # Assign IDs for time shapes: (TODO: check if it is really valid to add the zero to the time shape)
    oivrfun = [(o,i,cumsum([0; o.T])/seq.DEF["RadiofrequencyRasterTime"]) for (_,o,i) in boirfun if isa(o.T, Vector{<:Number})]
    [push!(wishape, (w, wishape[end][2]+1)) for (o,i,w) in oivrfun if all([!(length(w) == length(ws) && w ≈ ws) for (ws,_) in wishape])]
    # The important vector for rf
    oiawmpvthrfun = [[o, i, a, w, 0, 0, v, 0, 0] for (o,i,a,w,v) in oiawvrfun]
    for oiawmpvth in oiawmpvthrfun
        for (w,i) in wishape
            if length(w) == length(oiawmpvth[4])
                if w ≈ abs.(oiawmpvth[4])
                    oiawmpvth[5] = i
                #elseif w ≈ -angle.(oiawmpvth[4])/2π
                #    oiawmpvth[6] = i
                elseif all(w - -angle.(oiawmpvth[4])/2π .≈ w[1] - -angle.(oiawmpvth[4][1])/2π)
                    oiawmpvth[6] = i
                    a = -angle.(oiawmpvth[4][1]) - w[1]*2π
                    if a < 0
                        a += 2π
                    end
                    oiawmpvth[9] = a
                end
            elseif length(w) == length(oiawmpvth[7])
                if w ≈ oiawmpvth[7]
                    oiawmpvth[6] = i
                end
            end
        end
    end

    # For GR:
    ogron = [g for s in seq for g in s.GR[:,1] if sum(abs.(g.A)) > 0]
    oigrun = [(ogron[1], 1)]
    [push!(oigrun, (o, oigrun[end][2]+1)) for o in ogron if all([!(o ≈ ou) for (ou,iu) in oigrun])]
    # For GX:
    bogxon = [(b, s.GR[1,1]) for (b,s) in enumerate(seq) if KomaMRICore.is_Gx_on(s)]
    boigxon = [(b,o,iu) for (b,o) in bogxon for (ou,iu) in oigrun if o ≈ ou]
    boigx = [[b, s.GR[1,1], 0] for (b,s) in enumerate(seq)]
    [boi[3] = io for boi in boigx for (bo,oo,io) in boigxon if boi[1] == bo]
    # For GY:
    bogyon = [(b, s.GR[2,1]) for (b,s) in enumerate(seq) if KomaMRICore.is_Gy_on(s)]
    boigyon = [(b,o,iu) for (b,o) in bogyon for (ou,iu) in oigrun if o ≈ ou]
    boigy = [[b, s.GR[2,1], 0] for (b,s) in enumerate(seq)]
    [boi[3] = io for boi in boigy for (bo,oo,io) in boigyon if boi[1] == bo]
    # For GZ:
    bogzon = [(b, s.GR[3,1]) for (b,s) in enumerate(seq) if KomaMRICore.is_Gz_on(s)]
    boigzon = [(b,o,iu) for (b,o) in bogzon for (ou,iu) in oigrun if o ≈ ou]
    boigz = [[b, s.GR[3,1], 0] for (b,s) in enumerate(seq)]
    [boi[3] = io for boi in boigz for (bo,oo,io) in boigzon if boi[1] == bo]

    # For ADC
    boadcon = [(b, s.ADC[1]) for (b,s) in enumerate(seq) if is_ADC_on(s)]
    boiadcun = [(boadcon[1][1], boadcon[1][2], 1)]
    [push!(boiadcun, (b, o, boiadcun[end][3]+1)) for (b,o) in boadcon if all([!(o ≈ ou) for (bu,ou,iu) in boiadcun])]
    boiadcon = [(b,o,iu) for (b,o) in boadcon for (bu,ou,iu) in boiadcun if o ≈ ou]
    boiadc = [[b, s.ADC[1], 0] for (b,s) in enumerate(seq)]
    [boi[3] = io for boi in boiadc for (bo,oo,io) in boiadcon if boi[1] == bo]

    @printf(fid, "# Format of blocks:\n")
    @printf(fid, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
    @printf(fid, "[BLOCKS]\n")
    isEXT = 0
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
        idRF, idGX, idGY, idGZ, idADC = boirf[i][3], boigx[i][3], boigy[i][3], boigz[i][3], boiadc[i][3]
        Printf.format(fid, Printf.Format(idFormatStr * "%3d %3d %3d %3d %3d %2d %2d\n"), i, bdr, idRF, idGX, idGY, idGZ, idADC, isEXT)
    end
    @printf(fid, "\n")

    if !isempty(oiawmpvthrfun)
        @printf(fid, "# Format of RF events:\n")
        @printf(fid, "# id amplitude mag_id phase_id time_shape_id delay freq phase\n")
        @printf(fid, "# ..        Hz   ....     ....          ....    us   Hz   rad\n")
        @printf(fid, "[RF]\n")
        for (o,i,a,w,m,p,v,t,h) in oiawmpvthrfun
            id, amp, magid, phaseid, timeshapeid = i, a, m, p, t
            delay = round(o.delay / seq.DEF["RadiofrequencyRasterTime"]) * seq.DEF["RadiofrequencyRasterTime"] * 1e6
            freq, phase = o.Δf, h
            @printf(fid, "%d %12g %d %d %d %g %g %g\n", id, amp, magid, phaseid, timeshapeid, delay, freq, phase)
        end
        @printf(fid, "\n")
    end

    oigrgradun = [(o,i) for (o,i) in oigrun if length(o.A) != 1]
    oiawvgrgradun = [(o, i, magsign(o.A)*γ, o.A / magsign(o.A), cumsum([0; o.T])/seq.DEF["GradientRasterTime"]) for (o,i) in oigrgradun]
    # Assign IDs for RF shapes
    [push!(wishape, (w, wishape[end][2]+1)) for (o,i,a,w) in oiawvgrgradun if all([!(length(w) == length(ws) && w ≈ ws) for (ws,_) in wishape])]
    # Assign IDs for GRAD shapes
    oivrfun = [(o,i,cumsum([0; o.T])/seq.DEF["GradientRasterTime"]) for (o,i) in oigrgradun if isa(o.T, Vector{<:Number})]
    [push!(wishape, (w, wishape[end][2]+1)) for (_,_,w) in oivrfun if all([!(length(w) == length(ws) && w ≈ ws) for (ws,_) in wishape])]
    # The important grad vector
    oiawmvtgrgradun = [[o, i, a, w, 0, v, 0] for (o,i,a,w,v) in oiawvgrgradun]
    for oiawmvt in oiawmvtgrgradun
        for (w,i) in wishape
            if length(w) == length(oiawmvt[4]) && w ≈ oiawmvt[4]
                oiawmvt[5] = i
            elseif length(w) == length(oiawmvt[6]) && w ≈ oiawmvt[6]
                oiawmvt[7] = i
            end
        end
    end
    winshapewrite = [(length(compress(w)) == length(w) ? w : compress(w), i, length(w)) for (w,i) in wishape]
    oigrtrapun = [(o,i) for (o,i) in oigrun if length(o.A) == 1]

    if !isempty(oiawmvtgrgradun)
        @printf(fid, "# Format of arbitrary gradients:\n")
        @printf(fid, "#   time_shape_id of 0 means default timing (stepping with grad_raster starting at 1/2 of grad_raster)\n")
        @printf(fid, "# id amplitude amp_shape_id time_shape_id delay\n") # do we need delay ???
        @printf(fid, "# ..      Hz/m       ..         ..          us\n")
        @printf(fid, "[GRADIENTS]\n")
        for (o,i,a,w,m,v,t) in oiawmvtgrgradun
            id, amp, ampid, timeid = i, a, m, t
            delay =  round(o.delay * 1e6)
            @printf(fid, "%d %12g %d %d %d\n", id, amp, ampid, timeid, delay)
        end
        @printf(fid, "\n")
    end

    if !isempty(oigrtrapun)
        @printf(fid, "# Format of trapezoid gradients:\n")
        @printf(fid, "# id amplitude rise flat fall delay\n")
        @printf(fid, "# ..      Hz/m   us   us   us    us\n")
        @printf(fid, "[TRAP]\n")
        for (o,i) in oigrtrapun
            id, amp, rise, flat, fall, delay = i, o.A * γ, o.rise * 1e6, o.T * 1e6, o.fall * 1e6, o.delay * 1e6
            @printf(fid, "%2d %12g %3d %4d %3d %3d\n", id, amp, rise, flat, fall, delay)
        end
        @printf(fid, "\n")
    end

    if !isempty(boiadcun)
        @printf(fid, "# Format of ADC events:\n")
        @printf(fid, "# id num dwell delay freq phase\n")
        @printf(fid, "# ..  ..    ns    us   Hz   rad\n")
        @printf(fid, "[ADC]\n")
        for (_,o,i) in boiadcun
            id, num, dwell, delay, freq, phase = i, o.N, o.T * 1e9 / (o.N - 1), (o.delay - 0.5*o.T/(o.N - 1)) * 1e6, o.Δf, o.ϕ
            @printf(fid, "%d %d %.0f %.0f %g %g\n", id, num, dwell, delay, freq, phase)
        end
        @printf(fid, "\n")
    end

    if !isempty(winshapewrite)
        @printf(fid, "# Sequence Shapes\n")
        @printf(fid, "[SHAPES]\n\n")
        for (w,i,n) in winshapewrite
            @printf(fid, "shape_id %d\n", i)
            @printf(fid, "num_samples %d\n", n)
            [@printf(fid, "%.9g\n", wi) for wi in w]
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

# Check equality
Base.isapprox(s1::Sequence, s2::Sequence) = begin
    if length(s1) != length(s2)
        return false
    end
    return all([s1.ADC[i] ≈ s2.ADC[i] && s1.RF[i] ≈ s2.RF[i] && s1.GR[i] ≈ s2.GR[i] && s1.DUR[i] ≈ s2.DUR[i] for i in 1:length(s1)])
end

seq1 = read_seq(filename)
println("seq ≈ seq1")
println(seq ≈ seq1)
