module KomaMRI

#IMPORT PACKAGES
using Reexport
@reexport using KomaMRICore
@reexport using KomaMRIPlots
import KomaMRICore: update_blink_window_progress!

#GUI
using Blink, Interact, AssetRegistry
using MAT

#MRIReco
using MRIReco
@reexport using MRIReco: reconstruction

#Reconstruction
include("reconstruction/Recon.jl")

#GUI
include("ui/ExportMATFunctions.jl")
include("KomaUI.jl")

export KomaUI

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = Pkg.project().version

end
