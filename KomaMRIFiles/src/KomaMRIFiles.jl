module KomaMRIFiles

using KomaMRIBase
# .mrd .jld2
using FileIO
# For reading .seq
using Scanf
# For reading .h5/.mat phantoms
using HDF5, MAT
# For writing .seq
using Printf, MD5
# For writing .phantom
using InteractiveUtils

using Reexport
using MRIFiles
@reexport using MRIFiles: ISMRMRDFile
@reexport using FileIO: save

include("Sequence/ReadPulseq.jl")
include("Sequence/WritePulseq.jl")
include("Phantom/JEMRIS.jl")
include("Phantom/MRiLab.jl")
include("Phantom/Phantom.jl")

# Pulseq
export read_seq, write_seq
# Phantom
export read_phantom_jemris, read_phantom_MRiLab, read_phantom, write_phantom

# Package version: KomaMRIFiles.__VERSION__
using Pkg
__VERSION__ = VersionNumber(
    Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"]
)

end # module KomaMRIFiles
