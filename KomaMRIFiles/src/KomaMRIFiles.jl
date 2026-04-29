module KomaMRIFiles

using KomaMRIBase
import KomaMRIBase:
    AbstractSpinSpan,
    ArbitraryAction,
    BlockPulseRF,
    FrequencyModulatedRF,
    SimpleAction,
    TimeCurve,
    TimeShapedGrad,
    TimeShapedRF,
    TrapezoidalGrad,
    UniformlySampledGrad,
    extension_type_header,
    get_EXT_type_from_symbol,
    get_RF_use_from_char,
    get_char_from_RF_use,
    get_dims,
    _shape_times,
    get_scale,
    get_pulseq_format,
    get_symbol_from_EXT_type,
    is_on,
    sort_motions!
using FileIO, HDF5, MAT, InteractiveUtils, Printf # IO related
using SHA, MD5 # Pulseq signature verification
using Reexport
using MRIFiles
import MRIFiles: insertNode
@reexport using MRIFiles: ISMRMRDFile
@reexport using FileIO: save

include("Sequence/Pulseq.jl")
include("Phantom/JEMRIS.jl")
include("Phantom/MRiLab.jl")
include("Phantom/Phantom.jl")

export PulseqSequenceData, read_seq, read_seq_data, write_seq, write_seq_data    # Pulseq
export read_phantom_jemris, read_phantom_MRiLab, read_phantom, write_phantom     # Phantom

end # module KomaMRIFiles
