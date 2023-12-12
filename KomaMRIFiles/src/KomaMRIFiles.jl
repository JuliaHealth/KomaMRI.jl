module KomaMRIFiles

using KomaMRIBase
using Scanf, FileIO, HDF5, MAT    # IO related

using Reexport
using MRIFiles
import MRIFiles: insertNode
@reexport using MRIFiles: ISMRMRDFile
@reexport using FileIO: save

include("Sequence/Pulseq.jl")
include("Phantom/JEMRIS.jl")
include("Phantom/MRiLab.jl")

export read_seq                                 # Pulseq
export read_phantom_jemris, read_phantom_MRiLab # Phantom

# Package version: KomaMRIFiles.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaMRIFiles
