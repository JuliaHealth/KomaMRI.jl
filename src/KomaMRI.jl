module KomaMRI

#IMPORT PACKAGES
#Reconstruction
using MRIReco, MRIFiles
@reexport using MRIReco: RawAcquisitionData, AcquisitionData, reconstruction
@reexport using MRIFiles: ISMRMRDFile

using Reexport
@reexport using KomaMRICore

#GUI
using Blink, Interact, PlotlyJS, AssetRegistry
@reexport using PlotlyJS: savefig

#GUI
!Blink.AtomShell.isinstalled() && Blink.AtomShell.install()
include("KomaUI.jl")

export KomaUI

end
