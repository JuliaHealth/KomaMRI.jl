module KomaMRI

using Reexport
@reexport using KomaMRICore
@reexport using KomaMRIFiles
@reexport using KomaMRIPlots

# GUI
using Blink, Interact, AssetRegistry
using MAT

# Reconstruction
using FFTW: fftshift, ifftshift, fft, ifft
include("reconstruction/Recon.jl")

# MRIReco
using MRIReco
@reexport using MRIReco: reconstruction

#GUI
include("ui/ExportMATFunctions.jl")
include("ui/ExportUIFunctions.jl")
include("KomaUI.jl")

# Export the UI and the observables
export KomaUI
export sys_ui, seq_ui, obj_ui, raw_ui, img_ui

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
