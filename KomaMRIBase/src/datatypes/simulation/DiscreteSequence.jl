"""
    seqd = DiscreteSequence(Gx, Gy, Gz, B1, Δf, ADC, t, Δt)

A sampled version of a Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteSequence is the struct used for simulation.

# Arguments
- `Gx`: (`::AbstractVector{T<:Real}`, `[T/m]`) x-gradient vector
- `Gy`: (`::AbstractVector{T<:Real}`, `[T/m]`) y-gradient vector
- `Gz`: (`::AbstractVector{T<:Real}`, `[T/m]`) z-gradient vector
- `B1`: (`::AbstractVector{Complex{T<:Real}}`, `[T]`) RF amplitude vector
- `Δf`: (`::AbstractVector{T<:Real}`, `[Hz]`) RF carrier frequency displacement vector
- `ADC`: (`::AbstractVector{Bool}`) ADC sample vector
- `t`: (`::AbstractVector{T<:Real}`, `[s]`) time vector
- `Δt`: (`::AbstractVector{T<:Real}`, `[s]`) delta time vector

# Returns
- `seqd`: (`::DiscreteSequence`) DiscreteSequence struct
"""
struct DiscreteSequence{T<:Real}
    Gx::AbstractVector{T}
    Gy::AbstractVector{T}
    Gz::AbstractVector{T}
    B1::AbstractVector{Complex{T}}
    Δf::AbstractVector{T}
    ADC::AbstractVector{Bool}
    t::AbstractVector{T}
    Δt::AbstractVector{T}
end

Base.length(seq::DiscreteSequence) = length(seq.Δt)
Base.getindex(seq::DiscreteSequence, i::Integer) = begin
    DiscreteSequence(seq.Gx[i, :],
                     seq.Gy[i, :],
                     seq.Gz[i, :],
                     seq.B1[i, :],
                     seq.Δf[i, :],
                     seq.ADC[i, :],
                     seq.t[i, :],
                     seq.Δt[i, :])
end
Base.getindex(seq::DiscreteSequence, i::UnitRange) = begin
    DiscreteSequence(seq.Gx[i.start:i.stop+1],
                     seq.Gy[i.start:i.stop+1],
                     seq.Gz[i.start:i.stop+1],
                     seq.B1[i.start:i.stop+1],
                     seq.Δf[i.start:i.stop+1],
                     seq.ADC[i],
                     seq.t[i.start:i.stop+1],
                     seq.Δt[i])
end
Base.view(seq::DiscreteSequence, i::UnitRange) = begin
    @views DiscreteSequence(seq.Gx[i.start:i.stop+1],
                     seq.Gy[i.start:i.stop+1],
                     seq.Gz[i.start:i.stop+1],
                     seq.B1[i.start:i.stop+1],
                     seq.Δf[i.start:i.stop+1],
                     seq.ADC[i],
                     seq.t[i.start:i.stop+1],
                     seq.Δt[i])
end
Base.iterate(seq::DiscreteSequence) = (seq[1], 2)
Base.iterate(seq::DiscreteSequence, i) = (i <= length(seq)) ? (seq[i], i+1) : nothing

is_GR_on(seq::DiscreteSequence) =  sum(abs.([seq.Gx[1:end-1]; seq.Gy[1:end-1]; seq.Gz[1:end-1]])) != 0
is_RF_on(seq::DiscreteSequence) =  sum(abs.(seq.B1[1:end-1])) != 0
is_ADC_on(seq::DiscreteSequence) = sum(abs.(seq.ADC[1:end-1])) != 0
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
    B1, Δf     = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    tadc_set = Set(tadc)
    ADCflag = [tt in tadc_set for tt in t] # Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
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
