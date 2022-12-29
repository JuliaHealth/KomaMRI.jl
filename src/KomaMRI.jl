module KomaMRI

#IMPORT PACKAGES
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
