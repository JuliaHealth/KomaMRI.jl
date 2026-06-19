"""
    sys = Scanner(B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
        RF_ring_down_time, RF_dead_time, ADC_dead_time)

The Scanner struct. It contains hardware limitations of the MRI resonator. It is an input
for the simulation.

# Arguments
- `B0`: (`=1.5`, `[T]`) main magnetic field strength
- `B1`: (`=10e-6`, `[T]`) maximum RF amplitude
- `Gmax`: (`=60e-3`, `[T/m]`) maximum gradient amplitude
- `Smax`: (`=500.0`, `[mT/m/ms]`) gradient's maximum slew-rate
- `ADC_Δt`: (`=2e-6`, `[s]`) ADC raster time
- `DUR_Δt`: (`=1e-5`, `[s]`) block duration raster time
- `GR_Δt`: (`=1e-5`, `[s]`) gradient raster time
- `RF_Δt`: (`=1e-6`, `[s]`) RF raster time
- `RF_ring_down_time`: (`=20e-6`, `[s]`) RF ring down time
- `RF_dead_time`: (`=100e-6`, `[s]`) RF dead time
- `ADC_dead_time`: (`=10e-6`, `[s]`) ADC dead time

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
    B0::Float64 = 1.5
    B1::Float64 = 10e-6
    Gmax::Float64 = 60e-3
    Smax::Float64 = 500.0
    #Sampling
    ADC_Δt::Float64 = 2e-6
    DUR_Δt::Float64 = 1e-5
    GR_Δt::Float64 = 1e-5
    RF_Δt::Float64 = 1e-6
    #Secondary
    RF_ring_down_time::Float64 = 20e-6
    RF_dead_time::Float64 = 100e-6
    ADC_dead_time::Float64 = 10e-6

    function Scanner(
        B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
        RF_ring_down_time, RF_dead_time, ADC_dead_time,
    )
        return new(
            scanner_field_strength(B0), scanner_rf_amplitude(B1),
            scanner_gradient_amplitude(Gmax), scanner_slew_rate(Smax),
            scanner_time(ADC_Δt), scanner_time(DUR_Δt), scanner_time(GR_Δt),
            scanner_time(RF_Δt), scanner_time(RF_ring_down_time),
            scanner_time(RF_dead_time), scanner_time(ADC_dead_time),
        )
    end
end

scanner_field_strength(x) = Float64(x)
scanner_rf_amplitude(x) = Float64(x)
scanner_gradient_amplitude(x) = Float64(x)
scanner_slew_rate(x) = Float64(x)
scanner_time(x) = Float64(x)
