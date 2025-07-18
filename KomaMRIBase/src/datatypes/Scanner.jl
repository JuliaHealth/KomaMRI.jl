"""
    sys = Scanner(B0, B1, Gmax, Smax, ADC_Δt, seq_Δt, GR_Δt, RF_Δt,
        RF_ring_down_T, RF_dead_time_T, ADC_dead_time_T)

The `Scanner` struct contains the hardware limitations of an MRI resonator and serves as an
input for the simulation.

# Arguments
- `B0`: (`::Real`, `=1.5`, `[T]`) main magnetic field strength
- `B1`: (`::Real`, `=10e-6`, `[T]`) maximum RF amplitude
- `Gmax`: (`::Real`, `=60e-3`, `[T/m]`) maximum gradient amplitude
- `Smax`: (`::Real`, `=500`, `[mT/m/ms]`) gradient's maximum slew-rate
- `ADC_Δt`: (`::Real`, `=2e-6`, `[s]`) ADC raster time
- `seq_Δt`: (`::Real`, `=1e-5`, `[s]`) sequence-block raster time
- `GR_Δt`: (`::Real`, `=1e-5`, `[s]`) gradient raster time
- `RF_Δt`: (`::Real`, `=1e-6`, `[s]`) RF raster time
- `RF_ring_down_T`: (`::Real`, `=20e-6`, `[s]`) RF ring down time
- `RF_dead_time_T`: (`::Real`, `=100e-6`, `[s]`) RF dead time
- `ADC_dead_time_T`: (`::Real`, `=10e-6`, `[s]`) ADC dead time

# Returns
- `sys`: (`::Scanner`) Scanner struct

# Examples
```julia-repl
julia> sys = Scanner()

julia> sys.B0
```
"""
@with_kw mutable struct Scanner
    #Main
    B0::Real=1.5
    B1::Real=10e-6
    Gmax::Real=60e-3
    Smax::Real=500
    #Sampling
    ADC_Δt::Real=2e-6
    seq_Δt::Real=1e-5
    GR_Δt::Real=1e-5
    RF_Δt::Real=1e-6
    #Secondary
    RF_ring_down_T::Real=20e-6
    RF_dead_time_T::Real=100e-6
    ADC_dead_time_T::Real=10e-6
end
