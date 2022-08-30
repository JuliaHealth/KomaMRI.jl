"""
    scanner = Scanner(B0, B1, Gmax, Smax, ADC_Δt, seq_Δt, GR_Δt, RF_Δt,
        RF_ring_down_T, RF_dead_time_T, ADC_dead_time_T)

The Scanner struct.

# Arguments
- `B0`: (`::Real`, `=1.5`, `[T]`) the main magnetic field
- `B1`: (`::Real`, `=10e-6`, `[T]`) the maximum RF amplitude
- `Gmax`: (`::Real`, `=60e-3`, `[T/m]`) the maximum Gradient
- `Smax`: (`::Real`, `=500`, `[mT/m/ms]`) the maximum slew-rate
- `ADC_Δt`: (`::Real`, `=2e-6`, `[s]`) the ADC raster time
- `seq_Δt`: (`::Real`, `=1e-5`, `[s]`) the sequence-block raster time
- `GR_Δt`: (`::Real`, `=1e-5`, `[s]`) the gradient raster time
- `RF_Δt`: (`::Real`, `=1e-6`, `[s]`) the RF raster time
- `RF_ring_down_T`: (`::Real`, `=20e-6`, `[s]`) the RF ring down time
- `RF_dead_time_T`: (`::Real`, `=100e-6`, `[s]`) the RF dead time
- `ADC_dead_time_T`: (`::Real`, `=10e-6`, `[s]`) the ADC dead time

# Returns
- `scanner`: (`::Scanner`) the Scanner struct

# Examples
```julia-repl
julia> sys = Scanner()
Scanner
  B0: Float64 1.5
  B1: Float64 1.0e-5
  Gmax: Float64 0.06
  Smax: Int64 500
  ADC_Δt: Float64 2.0e-6
  seq_Δt: Float64 1.0e-5
  GR_Δt: Float64 1.0e-5
  RF_Δt: Float64 1.0e-6
  RF_ring_down_T: Float64 2.0e-5
  RF_dead_time_T: Float64 0.0001
  ADC_dead_time_T: Float64 1.0e-5
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

#Functions that check that Sequence satisfy hardware requirements:
#check_sys_req(seq::Sequence,sys::Scanner)
