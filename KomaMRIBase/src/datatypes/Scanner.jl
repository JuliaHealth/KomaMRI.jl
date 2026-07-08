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
- `gradients`: (`=LinearXYZ()`) gradient coil model
- `rf_tx`: (`=UniformTransmit()`) transmit coil model
- `rf_rx`: (`=UniformCoilSens()`) receive coil sensitivity model

# Returns
- `sys`: (`::Scanner`) Scanner struct

# Examples
```julia-repl
julia> sys = Scanner()

julia> sys.B0
```
"""
abstract type GradientCoils{T} end
struct LinearXYZ{T} <: GradientCoils{T} end

abstract type RFCoilsTx{T} end
struct UniformTransmit{T} <: RFCoilsTx{T} end

abstract type RFCoilsRx{T} end
struct UniformCoilSens{T} <: RFCoilsRx{T} end

Base.@kwdef struct TheoreticalCoilSensitivities{T} <: RFCoilsRx{T}
    ncoils::Int = 8
    radius::T = T(0.20)
    L::T = T(0.30)
end

export HardwareLimits, Scanner, GradientCoils, LinearXYZ, RFCoilsTx, UniformTransmit,
    RFCoilsRx, UniformCoilSens, TheoreticalCoilSensitivities, get_n_coils

get_n_coils(::RFCoilsRx) = 1
get_n_coils(coils::TheoreticalCoilSensitivities) = coils.ncoils

coil_sensitivities(
    x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, ::UniformCoilSens
) where {T<:Real} = ones(Complex{T}, 1, length(x))

@with_kw mutable struct HardwareLimits
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

    function HardwareLimits(
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

mutable struct Scanner{
    Limits <: HardwareLimits,
    RFTx <: RFCoilsTx,
    RFRx <: RFCoilsRx,
    GRCoils <: GradientCoils,
}
    limits::Limits
    rf_tx::RFTx
    rf_rx::RFRx
    gradients::GRCoils
end

function Scanner(
    B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
    RF_ring_down_time, RF_dead_time, ADC_dead_time;
    gradients=LinearXYZ{Float64}(),
    rf_tx=UniformTransmit{Float64}(),
    rf_rx=UniformCoilSens{Float64}(),
    rf_coils=nothing,
)
    isnothing(rf_coils) || (rf_rx = rf_coils)
    limits = HardwareLimits(
        B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
        RF_ring_down_time, RF_dead_time, ADC_dead_time,
    )
    return Scanner{typeof(limits),typeof(rf_tx),typeof(rf_rx),typeof(gradients)}(
        limits,
        rf_tx,
        rf_rx,
        gradients,
    )
end

function Scanner(;
    B0=1.5,
    B1=10e-6,
    Gmax=60e-3,
    Smax=500.0,
    ADC_Δt=2e-6,
    DUR_Δt=1e-5,
    GR_Δt=1e-5,
    RF_Δt=1e-6,
    RF_ring_down_time=20e-6,
    RF_dead_time=100e-6,
    ADC_dead_time=10e-6,
    gradients=LinearXYZ{Float64}(),
    rf_tx=UniformTransmit{Float64}(),
    rf_rx=UniformCoilSens{Float64}(),
    rf_coils=nothing,
)
    isnothing(rf_coils) || (rf_rx = rf_coils)
    return Scanner(
        B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
        RF_ring_down_time, RF_dead_time, ADC_dead_time;
        gradients, rf_tx, rf_rx,
    )
end

function coil_sensitivities(
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T};
    ncoils=8,
    radius=T(0.20),
    L=T(0.30)
) where {T<:Real}
    nspins = length(x)
    sens = Matrix{Complex{T}}(undef, ncoils, nspins)

    ϵ = eps(T)

    for n in 1:ncoils
        ϕn = T(2π) * T(n - 1) / T(ncoils)

        xn = radius * cos(ϕn)
        yn = radius * sin(ϕn)

        for j in 1:nspins
            rn = sqrt((x[j] - xn)^2 + (y[j] - yn)^2 + ϵ)

            Bϕ =
                (1 / rn) *
                (
                    (L - z[j]) / sqrt(rn^2 + (L - z[j])^2) +
                    (L + z[j]) / sqrt(rn^2 + (L + z[j])^2)
                )

            Bx = -Bϕ * (y[j] - yn) / rn
            By = Bϕ * (x[j] - xn) / rn

            sens[n, j] = Complex{T}(0.5) * (Bx - im * By)
        end

        sens[n, :] ./= maximum(abs.(sens[n, :])) + ϵ
    end
    return sens
end

function coil_sensitivities(
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    rf_rx::TheoreticalCoilSensitivities,
) where {T<:Real}
    return coil_sensitivities(
        x, y, z;
        ncoils=rf_rx.ncoils,
        radius=T(rf_rx.radius),
        L=T(rf_rx.L),
    )
end

function Base.getproperty(sys::Scanner, key::Symbol)
    if key === :rf_sens
        rf_rx = getfield(sys, :rf_rx)
        return (x, y, z) -> coil_sensitivities(x, y, z, rf_rx)
    elseif key === :rf_coils
        return getfield(sys, :rf_rx)
    elseif key in fieldnames(HardwareLimits)
        return getfield(getfield(sys, :limits), key)
    else
        return getfield(sys, key)
    end
end

function Base.setproperty!(sys::Scanner, key::Symbol, value)
    if key === :rf_coils
        setfield!(sys, :rf_rx, value)
    elseif key === :B0
        setfield!(getfield(sys, :limits), :B0, scanner_field_strength(value))
    elseif key === :B1
        setfield!(getfield(sys, :limits), :B1, scanner_rf_amplitude(value))
    elseif key === :Gmax
        setfield!(getfield(sys, :limits), :Gmax, scanner_gradient_amplitude(value))
    elseif key === :Smax
        setfield!(getfield(sys, :limits), :Smax, scanner_slew_rate(value))
    elseif key === :ADC_Δt
        setfield!(getfield(sys, :limits), :ADC_Δt, scanner_time(value))
    elseif key === :DUR_Δt
        setfield!(getfield(sys, :limits), :DUR_Δt, scanner_time(value))
    elseif key === :GR_Δt
        setfield!(getfield(sys, :limits), :GR_Δt, scanner_time(value))
    elseif key === :RF_Δt
        setfield!(getfield(sys, :limits), :RF_Δt, scanner_time(value))
    elseif key === :RF_ring_down_time
        setfield!(getfield(sys, :limits), :RF_ring_down_time, scanner_time(value))
    elseif key === :RF_dead_time
        setfield!(getfield(sys, :limits), :RF_dead_time, scanner_time(value))
    elseif key === :ADC_dead_time
        setfield!(getfield(sys, :limits), :ADC_dead_time, scanner_time(value))
    else
        setfield!(sys, key, value)
    end
    return value
end

scanner_field_strength(x) = Float64(x)
scanner_rf_amplitude(x) = Float64(x)
scanner_gradient_amplitude(x) = Float64(x)
scanner_slew_rate(x) = Float64(x)
scanner_time(x) = Float64(x)
