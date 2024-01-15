module KomaMRIFiles

using KomaMRIBase
using Scanf, FileIO, HDF5, MAT    # IO related
using Printf, MD5                 # For writing .seq files from Sequence struct

using Reexport
using MRIFiles
import MRIFiles: insertNode
@reexport using MRIFiles: ISMRMRDFile
@reexport using FileIO: save

include("Sequence/ReadPulseq.jl")
include("Sequence/WritePulseq.jl")
include("Phantom/JEMRIS.jl")
include("Phantom/MRiLab.jl")

export read_seq, write_seq                      # Pulseq
export read_phantom_jemris, read_phantom_MRiLab # Phantom

# Package version: KomaMRIFiles.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaMRIFiles
