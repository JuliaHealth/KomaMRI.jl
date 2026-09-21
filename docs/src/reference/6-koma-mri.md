# KomaMRI

```@meta
CurrentModule = KomaMRI
```

## User Interface

```@docs
KomaUI
load_file!
```

### Actions

```@docs
click!
```

### Window data

Each window owns its data and options:

```julia
w.seq[] = seq        # Sequence
w.obj[] = obj        # Phantom
w.sys[] = sys        # Scanner
w.physio[] = physio  # Physiological signal
w.raw[] = raw        # Raw data
w.img[] = img        # Image
w.sim_params[] = sim_params  # Dict{String,Any}
w.rec_params[] = rec_params  # Dict{Symbol,Any}
```

Assignments refresh that window's view without changing filenames or reload paths.
In-place edits require `notify(observable)` to refresh the view. Simulation,
reconstruction, and export read the current values.

Sequence assignment resets the physiological signal to its default. Scanner updates
show coil sensitivities when the receiver changes, otherwise hardware limits when
limits change. With neither changed, the current Scanner view refreshes, or coil
sensitivities open if another section is displayed.

`Scanner` is immutable; receiver replacement preserves the other components:

```julia
sys = w.sys[]
w.sys[] = Scanner(; limits=sys.limits, gradient=sys.gradient,
    transmitter=sys.transmitter, receiver=BirdcageCoilSens())
```
