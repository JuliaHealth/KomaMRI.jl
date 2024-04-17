abstract type MRISequenceEvent end #get all available types by using subtypes(KomaMRIBase.MRISequenceEvent)

function Base.:≈(ev1::T, ev2::T; kwargs...) where {T<:MRISequenceEvent}
    if any(length(getfield(ev1, k)) != length(getfield(ev2, k)) for k in fieldnames(T))
        return false
    end
    return all(≈(getfield(ev1, k), getfield(ev2, k); kwargs...) for k in fieldnames(T))
end

function Base.:(==)(ev1::T, ev2::T) where {T<:MRISequenceEvent}
    if any(length(getfield(ev1, k)) != length(getfield(ev2, k)) for k in fieldnames(T))
        return false
    end
    return all(getfield(ev1, k) == getfield(ev2, k) for k in fieldnames(T))
end

Base.size(ev::T, i) where {T<:MRISequenceEvent} = 1

# dur(ev::AbstractArray{T}) where {T<:MRISequenceEvent} = maximum(dur.(ev), dims=1)[:]
