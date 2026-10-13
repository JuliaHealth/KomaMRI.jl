# KomaMRIPlots

```@meta
CurrentModule = KomaMRIPlots
```

Plots display directly in Pluto without an extra wrapper or plotting dependency.
Reruns preserve valid zoom, camera and legend selections.
Nonadaptive plots keep their `PlotlyBase.Plot` return type and Plotly mutation/export APIs.

Adaptive plots require live Julia in REPL, Pluto, Jupyter, Quarto, and KomaUI.
Nonadaptive figures work offline unless `adaptive_navigator=true` adds a live
sequence panel to k-space. For Quarto live views, keep the Julia worker running
with `execute.daemon`.

## Plotting `Phantom`

```@docs
plot_phantom_map
SpatialPlot
plot_phantom
plot_phantom_slice
```

## Plotting receive sensitivities

```@docs
plot_coil_sens
get_coil_sens_fov
```

## Plotting `Sequence`

Sequence-like plots and `plot_signal` return a regular `PlotlyBase.Plot` by default.
Set `adaptive=true` to fetch detailed samples on zoom or pan. This live mode requires
the Julia session to remain running; KomaUI enables it automatically.

```julia
plot_seq(seq)                  # Regular Plotly figure
plot_seq(seq; adaptive=true)   # Live adaptive viewer
```

```@docs
TimePlot
plot_seq
plot_kspace
plot_M0
plot_M1
plot_M2
plot_eddy_currents
plot_slew_rate
plot_seqd
```

## Plotting `RawAcquisitionData`

```@docs
plot_signal
```

## Plotting images

```@docs
plot_image
```

## Others

```@docs
plot_dict
savefig
```
