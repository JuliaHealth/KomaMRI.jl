module KomaMRI

# IMPORT PACKAGES
using Reexport
@reexport using KomaMRICore
@reexport using KomaMRIFiles
@reexport using KomaMRIPlots

using Bonito
using Artifacts
import Electron
import MsgPack
import PlotlyBase
using Observables: Observable, ObserverFunction, off
using MAT

# Reconstruction
using FFTW: fftshift, ifftshift, fft, ifft
include("reconstruction/Recon.jl")

# MRIReco
using MRIReco
@reexport using MRIReco: reconstruction

# UI
include("ui/ExportMATFunctions.jl")
include("ui/UIDefaults.jl")
include("KomaUI.jl")
include("KomaCLI.jl")
@static if VERSION >= v"1.12"
    include("KomaApp.jl")
end

# Export the UI and the observables
export KomaUI, KomaWindow
export sys_ui, seq_ui, obj_ui, raw_ui, img_ui

end
