# KomaMRI

```@meta
CurrentModule = KomaMRI
```

## User Interface

```@docs
KomaUI
KomaWorkspace
add_window!
load_file!
simulate!
reconstruct!
export_results
ui_status
```

### Actions

```@docs
click!
ui_actions
```

### Window data

Each window owns its data and parameters:

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

The `:recon_policy` entry in `w.rec_params[]` controls label handling and coil
combination, e.g. `(COIL=:rss,)`. It is displayed with the reconstruction parameters
and applied by `reconstruct!`.

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

## Reconstruction

```@docs
reconstruct_with_labels
reconstruction_policy
DEFAULT_RECON_POLICY
```
