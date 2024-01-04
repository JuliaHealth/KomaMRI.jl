module KomaMRIPlots

using KomaMRIBase
using MAT, Interpolations, PlotlyJS

include("ui/PlotBackends.jl")
include("ui/DisplayFunctions.jl")

using Reexport
@reexport using PlotlyJS: savefig

export plot_seq, plot_M0, plot_M1, plot_M2, plot_eddy_currents, plot_seqd,
        plot_slew_rate, plot_kspace, plot_phantom_map, plot_signal, plot_image, plot_dict

#Package version, KomaMRIPlots.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end
