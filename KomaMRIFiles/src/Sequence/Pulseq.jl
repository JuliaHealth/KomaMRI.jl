include("pulseq/PulseqEvents.jl")

const QUANT_TOL = 1e-12
const DEFAULT_RASTER = PulseqRaster(Sequence(), Scanner())

include("pulseq/Signature.jl")

include("pulseq/WritePulseq.jl")
include("pulseq/ReadPulseq.jl")