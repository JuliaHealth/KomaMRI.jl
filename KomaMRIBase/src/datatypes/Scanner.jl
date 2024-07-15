# Hardware limits
@with_kw mutable struct HardwareLimits{T}
    B0::T = 1.5
    B1::T = 10e-6
    Gmax::T = 60e-3
    Smax::T = 500.0
    ADC_Δt::T = 2e-6
    seq_Δt::T = 1e-5
    GR_Δt::T = 1e-5
    RF_Δt::T = 1e-6
    RF_ring_down_T::T = 20e-6
    RF_dead_time_T::T = 100e-6
    ADC_dead_time_T::T = 10e-6
end

# Gradients
abstract type Gradients{T} end
@with_kw mutable struct LinearXYZGradients{T} <: Gradients{T} 
    Gx::AbstractVector{T} = zeros(T, 1)
    Gy::AbstractVector{T} = zeros(T, 1)
    Gz::AbstractVector{T} = zeros(T, 1)
end

# RF coils
abstract type RFCoils{T} end
@with_kw mutable struct UniformRFCoils{T} <: RFCoils{T} 
    coil_sens::AbstractMatrix{Complex{T}} = complex(ones(Complex{T}, 1, 1))
end

struct ArbitraryRFCoils{T} <: RFCoils{T}
    x::AbstractVector{T} 
    y::AbstractVector{T} 
    z::AbstractVector{T} 
    coil_sens::AbstractMatrix{Complex{T}}  
    B1⁺::AbstractMatrix{Complex{T}} 
end

struct RFCoilsSensDefinedAtPhantomPositions{T} <: RFCoils{T}
    coil_sens::AbstractMatrix{Complex{T}}
end

"""
    sys = Scanner(B0, B1, Gmax, Smax, ADC_Δt, seq_Δt, GR_Δt, RF_Δt,
        RF_ring_down_T, RF_dead_time_T, ADC_dead_time_T)

The Scanner struct. It contains hardware limitations of the MRI resonator. It is an input
for the simulation.

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
@with_kw mutable struct Scanner{T}
    limits::HardwareLimits{T} = HardwareLimits{Float64}()
    gradients::Gradients{T} = LinearXYZGradients{Float64}()
    rf_coils::RFCoils{T} = UniformRFCoils{Float64}()
end

function Base.view(sys::Scanner, p)
    return Scanner(limits=sys.limits, gradients=sys.gradients, rf_coils=@view(sys.rf_coils[p]))
end

function Base.view(rf_coils::RFCoilsSensDefinedAtPhantomPositions, p)
    return RFCoilsSensDefinedAtPhantomPositions(@view rf_coils.coil_sens[p,:])
end

function Base.view(rf_coils::UniformRFCoils, p)
    return rf_coils
end

function acquire_signal!(sig, rf_coils::UniformRFCoils, Mxy)
    sig .= transpose(sum(Mxy; dims=1))
    return nothing
end


function acquire_signal!(sig, rf_coils::RFCoilsSensDefinedAtPhantomPositions, Mxy)
    for i in 1:size(rf_coils.coil_sens, 2)
        sig[:, i] .= transpose(sum(rf_coils.coil_sens[:, i] .* Mxy; dims=1))
    end
    return nothing
end

export ArbitraryRFCoils, RFCoilsSensDefinedAtPhantomPositions, UniformRFCoils, acquire_signal!, HardwareLimits, LinearXYZGradients
