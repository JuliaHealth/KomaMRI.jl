# KomaMRIPlots

```@meta
CurrentModule = KomaMRIPlots
```

## Plotting `Phantom`

```@docs
plot_phantom_map
SpatialPlot
plot_phantom
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
