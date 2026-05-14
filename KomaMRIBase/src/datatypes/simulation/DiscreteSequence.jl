"""
    seqd = DiscreteSequence(Gx, Gy, Gz, B1, Δf, ψ, ADC, t, Δt)

A sampled version of a Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteSequence is the struct used for simulation.

# Arguments
- `Gx`: (`::AbstractVector{T<:Real}`, `[T/m]`) x-gradient vector
- `Gy`: (`::AbstractVector{T<:Real}`, `[T/m]`) y-gradient vector
- `Gz`: (`::AbstractVector{T<:Real}`, `[T/m]`) z-gradient vector
- `B1`: (`::AbstractVector{Complex{T<:Real}}`, `[T]`) RF amplitude vector
- `Δf`: (`::AbstractVector{T<:Real}`, `[Hz]`) RF carrier frequency displacement vector
- `ψ`: (`::AbstractVector{T<:Real}`, `[rad]`) RF rotating-frame phase vector
- `ADC`: (`::AbstractVector{Bool}`) ADC sample vector
- `t`: (`::AbstractVector{T<:Real}`, `[s]`) time vector
- `Δt`: (`::AbstractVector{T<:Real}`, `[s]`) delta time vector

# Returns
- `seqd`: (`::DiscreteSequence`) DiscreteSequence struct
"""
struct DiscreteSequence{T<:Real,GxType,GyType,GzType,B1Type,ΔfType,ψType,ADCType,tType,ΔtType}
    Gx::GxType
    Gy::GyType
    Gz::GzType
    B1::B1Type
    Δf::ΔfType
    ψ::ψType
    ADC::ADCType
    t::tType
    Δt::ΔtType
end

function DiscreteSequence(
    Gx::AbstractVector{T},
    Gy::AbstractVector{T},
    Gz::AbstractVector{T},
    B1::AbstractVector{Complex{T}},
    Δf::AbstractVector{T},
    ψ::AbstractVector{T},
    ADC::AbstractVector{Bool},
    t::AbstractVector{T},
    Δt::AbstractVector{T},
) where {T<:Real}
    return DiscreteSequence{
        T,
        typeof(Gx),
        typeof(Gy),
        typeof(Gz),
        typeof(B1),
        typeof(Δf),
        typeof(ψ),
        typeof(ADC),
        typeof(t),
        typeof(Δt),
    }(Gx, Gy, Gz, B1, Δf, ψ, ADC, t, Δt)
end

function DiscreteSequence(Gx, Gy, Gz, B1, Δf, ψ, ADC, t, Δt)
    T = promote_type(
        _real_storage_eltype(typeof(Gx)),
        _real_storage_eltype(typeof(Gy)),
        _real_storage_eltype(typeof(Gz)),
        _real_storage_eltype(typeof(B1)),
        _real_storage_eltype(typeof(Δf)),
        _real_storage_eltype(typeof(ψ)),
        _real_storage_eltype(typeof(t)),
        _real_storage_eltype(typeof(Δt)),
    )
    return DiscreteSequence{
        T,
        typeof(Gx),
        typeof(Gy),
        typeof(Gz),
        typeof(B1),
        typeof(Δf),
        typeof(ψ),
        typeof(ADC),
        typeof(t),
        typeof(Δt),
    }(Gx, Gy, Gz, B1, Δf, ψ, ADC, t, Δt)
end

function _storage_eltype(::Type{A}) where {A}
    params = Base.unwrap_unionall(A).parameters
    return !isempty(params) && params[1] isa Type ? params[1] : eltype(A)
end
_real_storage_eltype(::Type{A}) where {A} = _real_eltype(_storage_eltype(A))
_real_eltype(::Type{Complex{T}}) where {T} = T
_real_eltype(::Type{T}) where {T<:Real} = T

Base.length(seq::DiscreteSequence) = length(seq.Δt)
Base.getindex(seq::DiscreteSequence, i::Integer) = begin
    DiscreteSequence(seq.Gx[i, :],
                     seq.Gy[i, :],
                     seq.Gz[i, :],
                     seq.B1[i, :],
                     seq.Δf[i, :],
                     seq.ψ[i, :],
                     seq.ADC[i, :],
                     seq.t[i, :],
                     seq.Δt[i, :])
end
Base.getindex(seq::DiscreteSequence, i::UnitRange) = begin
    DiscreteSequence(seq.Gx[i],
                     seq.Gy[i],
                     seq.Gz[i],
                     seq.B1[i],
                     seq.Δf[i],
                     seq.ψ[i],
                     seq.ADC[i],
                     seq.t[i],
                     seq.Δt[i.start:i.stop-1])
end
Base.view(seq::DiscreteSequence, i::UnitRange) = begin
    @views DiscreteSequence(seq.Gx[i],
                     seq.Gy[i],
                     seq.Gz[i],
                     seq.B1[i],
                     seq.Δf[i],
                     seq.ψ[i],
                     seq.ADC[i],
                     seq.t[i],
                     seq.Δt[i.start:i.stop-1])
end
Base.iterate(seq::DiscreteSequence) = (seq[1], 2)
Base.iterate(seq::DiscreteSequence, i) = (i <= length(seq)) ? (seq[i], i+1) : nothing

is_GR_on(seq::DiscreteSequence) =  sum(abs.([seq.Gx; seq.Gy; seq.Gz])) != 0
is_RF_on(seq::DiscreteSequence) =  sum(abs.(seq.B1)) != 0
is_ADC_on(seq::DiscreteSequence) = sum(abs.(seq.ADC)) != 0
is_GR_off(seq::DiscreteSequence) =  !is_GR_on(seq)
is_RF_off(seq::DiscreteSequence) =  !is_RF_on(seq)
is_ADC_off(seq::DiscreteSequence) = !is_ADC_on(seq)

"""
    seqd = discretize(seq::Sequence; sampling_params=default_sampling_params())

This function returns a sampled Sequence struct with RF and gradient time refinements
based on simulation parameters.

# Arguments
- `seq`: (`::Sequence`) sequence

# Keywords
- `sampling_params`: (`::Dict{String, Any}`, `=default_sampling_params()`) sampling
    parameter dictionary

# Returns
- `seqd`: (`::DiscreteSequence`) DiscreteSequence struct
"""
function discretize(seq::Sequence; sampling_params=default_sampling_params(), motion=NoMotion())
    t, Δt      = get_variable_times(seq; Δt=sampling_params["Δt"], Δt_rf=sampling_params["Δt_rf"], motion=motion)
    B1, Δf, ψ  = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    tadc_set = Set(tadc)
    ADCflag = [tt in tadc_set for tt in t] # Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ψ, ADCflag, t, Δt)
    return seqd
end

"""
Returns a dictionary with default simulation parameters.
"""
function default_sampling_params(sampling_params=Dict{String,Any}())
    get!(sampling_params, "Δt", 1e-3)
    get!(sampling_params, "Δt_rf", 5e-5)
    return sampling_params
end
