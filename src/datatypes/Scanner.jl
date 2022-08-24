"""
    Scanner(args...)

The Scanner object.

# Arguments
- `B0::Real=1.5`: main magnetic field [T]
- `B1::Real=10e-6`: max RF amplitude [T]
- `Gmax::Real=60e-3`: max Gradient [T/m]
- `Smax::Real=500`: max slew-rate [mT/m/ms] or [T/m/s]
- `ADC_Δt::Real=2e-6`: ADC raster time [s]
- `seq_Δt::Real=1e-5`: Sequence-block raster time [s]
- `GR_Δt::Real=1e-5`: Gradient raster time [s]
- `RF_Δt::Real=1e-6`: RF raster time [s]
- `RF_ring_down_T::Real=20e-6`: RF ring down time [s]
- `RF_dead_time_T::Real=100e-6`: RF dead tim [s]
- `ADC_dead_time_T::Real=10e-6`: ADC dead time [s]

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
    B0::Real=1.5      # Main magnetic field [T]
    B1::Real=10e-6    # Max RF amplitude [T]
    Gmax::Real=60e-3  # Max Gradient [T/m]
    Smax::Real=500    # Max Slew-rate [mT/m/ms] or [T/m/s]
    #Sampling
    ADC_Δt::Real=2e-6 # ADC raster time [s]
    seq_Δt::Real=1e-5 # Seq-block raster time [s]
    GR_Δt::Real=1e-5  # GR raster time [s]
    RF_Δt::Real=1e-6  # RF raster time [s]
    #Secondary
    RF_ring_down_T::Real=20e-6  # RF ring down time [s]
    RF_dead_time_T::Real=100e-6 # RF dead tim [s]
    ADC_dead_time_T::Real=10e-6 # ADC dead time [s]
end

#Functions that check that Sequence satisfy hardware requirements:
#check_sys_req(seq::Sequence,sys::Scanner)
