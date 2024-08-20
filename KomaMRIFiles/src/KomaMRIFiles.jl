module KomaMRIFiles

using KomaMRIBase
using Scanf, FileIO, HDF5, MAT, InteractiveUtils # IO related

using Reexport
using MRIFiles
import MRIFiles: insertNode
@reexport using MRIFiles: ISMRMRDFile
@reexport using FileIO: save

include("Sequence/Pulseq.jl")
include("Phantom/JEMRIS.jl")
include("Phantom/MRiLab.jl")
include("Phantom/Phantom.jl")

export read_seq                                                                     # Pulseq
export read_phantom_jemris, read_phantom_MRiLab, read_phantom, write_phantom        # Phantom

# Package version: KomaMRIFiles.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaMRIFiles
