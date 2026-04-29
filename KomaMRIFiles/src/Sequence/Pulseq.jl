include("pulseq/PulseqEvents.jl")
include("pulseq/Signature.jl")
include("pulseq/ReadPulseq.jl")
include("pulseq/WritePulseq.jl")

# Mirrors `teps = 1e-12` in MATLAB Pulseq
# `matlab/+mr/@Sequence/Sequence.m` (v1.5.1, line 2198).
const PULSEQ_TIME_TOL = 1e-12
const DEFAULT_RASTER = PulseqRaster(KomaMRIBase.Sequence())

"""
compress_shape Compress a gradient or pulse shape.
   s=compress_shape(w) Compress the waveform using a run-length compression
   scheme on the derivative. This strategy encodes constant and linear
   waveforms with very few samples. Returns:
     num_samples - the number of samples in the uncompressed waveform
     data - containing the compressed waveform

   See also decompress_shape
"""
# Mirrors `quant_fac = 1e-7` in MATLAB Pulseq
# `matlab/+mr/compressShape.m` (v1.5.1, line 33).
const PULSEQ_SHAPE_QUANTIZATION = 1e-7
# Mirrors `v(abs(v)<1e-10)=0` in MATLAB Pulseq
# `matlab/+mr/compressShape.m` (v1.5.1, line 52).
const PULSEQ_SHAPE_ZERO_TOL = 1e-10

function compress_shape(w; forceCompression=false)
    num_samples = length(w)
    if !forceCompression && num_samples <= 4
        data = copy(w)
    else
        quant_fac = PULSEQ_SHAPE_QUANTIZATION
        ws = w ./ quant_fac
        datq = round.([ws[1]; diff(ws[:])])
        qerr = ws[:] .- cumsum(datq)
        qcor = [0; diff(round.(qerr))]
        datd = datq .+ qcor
        maskChanges = [true; diff(datd) .!= 0]
        vals = datd[maskChanges] .* quant_fac
        k = findall([maskChanges; true])
        n = diff(k)
        nExtra = n .- 2
        v = Float64[]
        sizehint!(v, length(vals) + 2 * count(x -> x >= 0, nExtra))
        for i in eachindex(vals)
            push!(v, vals[i])
            if nExtra[i] >= 0
                push!(v, vals[i])
                push!(v, Float64(nExtra[i]))
            end
        end
        v[abs.(v) .<= PULSEQ_SHAPE_ZERO_TOL] .= 0
        data = forceCompression || num_samples > length(v) ? v : copy(w)
    end
    return num_samples, data
end

"""
decompress_shape Decompress a gradient or pulse shape.
   w=decompress_shape(shape) Decompress the shape compressed with a run-length
   compression scheme on the derivative. Returns:
     num_samples - the number of samples in the uncompressed waveform
     data - containing the compressed waveform

   See also compress_shape
"""
function decompress_shape(num_samples, data; forceDecompression = false)
    dataPack = data
    dataPackLen = length(dataPack)
    numSamples = num_samples
    if !forceDecompression && numSamples == dataPackLen
        w = dataPack
    else
        w = zeros(numSamples)
        dataPackDiff = dataPack[2:end] .- dataPack[1:end-1]
        dataPackMarkers = findall(dataPackDiff .== 0)
        countPack = 1
        countUnpack = 1
        for i in eachindex(dataPackMarkers)
            nextPack = dataPackMarkers[i]
            currUnpackSamples = nextPack - countPack
            if currUnpackSamples < 0
                continue
            elseif currUnpackSamples > 0
                w[countUnpack:(countUnpack + currUnpackSamples - 1)] .= dataPack[countPack:(nextPack - 1)]
                countPack += currUnpackSamples
                countUnpack += currUnpackSamples
            end
            rep = floor(Int, dataPack[countPack + 2] + 2)
            w[countUnpack:(countUnpack + rep - 1)] .= dataPack[countPack]
            countPack += 3
            countUnpack += rep
        end
        if countPack <= dataPackLen
            @assert dataPackLen - countPack == numSamples - countUnpack "Unsuccessful unpacking of samples"
            w[countUnpack:end] .= dataPack[countPack:end]
        end
        w = cumsum(w)
    end
    return w
end
