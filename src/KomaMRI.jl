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
const global CONT_INDEX = "index"
const global CONT_SEQUENCE = "sequence"
const global CONT_KSPACE = "kspace"
const global CONT_M0 = "m0"
const global CONT_M1 = "m1"
const global CONT_M2 = "m2"
const global CONT_PHANTOM = "phantom"
const global CONT_SCANNER_PARAMS = "scanneparams"
const global CONT_SIM_PARAMS = "simparams"
const global CONT_REC_PARAMS = "recparams"
const global CONT_ABSI = "absi"
const global CONT_ANGI = "angi"
const global CONT_ABSK = "absk"
const global CONT_SIG = "sig"
include("KomaUI.jl")
export CONT_INDEX, CONT_SEQUENCE, CONT_KSPACE, CONT_M0, CONT_M1, CONT_M2, CONT_PHANTOM,
    CONT_SCANNER_PARAMS, CONT_SIM_PARAMS, CONT_REC_PARAMS, CONT_ABSI, CONT_ANGI, CONT_ABSK, CONT_SIG

export KomaUI

end
