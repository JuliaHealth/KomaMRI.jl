module KomaMRI

#IMPORT PACKAGES
using Reexport
@reexport using KomaMRICore
import KomaMRICore: update_blink_window_progress!

#GUI
using Blink, Interact, PlotlyJS, AssetRegistry
@reexport using PlotlyJS: savefig

#GUI
!Blink.AtomShell.isinstalled() && Blink.AtomShell.install()
include("KomaUI.jl")

export KomaUI

end
