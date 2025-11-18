function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteSequence{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    @. Bz = x * seq.Gx[seq_idx] + y * seq.Gy[seq_idx] + z * seq.Gz[seq_idx]
end

function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteHigherOrderSequence{T, -1},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    @. Bz = x * seq.G[1, seq_idx] + y * seq.G[2, seq_idx] + z * seq.G[3, seq_idx]
    return nothing
end

function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteHigherOrderSequence{T, 0},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    Bz .= seq.G[1, seq_idx]
    return nothings
end

function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteHigherOrderSequence{T, 1},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    @. Bz = seq.G[1, seq_idx] + x * seq.G[2, seq_idx] + y * seq.G[3, seq_idx] + z * seq.G[4, seq_idx]
    return nothing
end

function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteHigherOrderSequence{T, 2},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    Bz .= seq.G[1, seq_idx]
    Bz .+= seq.G[2, seq_idx] .* x
    Bz .+= seq.G[3, seq_idx] .* y
    Bz .+= seq.G[4, seq_idx] .* z
    Bz .+= seq.G[5, seq_idx] .* x .* y # XY
    Bz .+= seq.G[6, seq_idx] .* y .* z # YZ
    Bz .+= seq.G[7, seq_idx] .* (2.0.*z.^2 - x.^2 - y.^2) # Z2
    Bz .+= seq.G[8, seq_idx] .* x .* z # XZ
    Bz .+= seq.G[9, seq_idx] .* (x.^2 - y.^2) # X2-Y2
    return nothing
end

function get_Bz_field!(
    Bz::AbstractVector{T},
    seq::DiscreteHigherOrderSequence,
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    seq_idx::Integer
) where {T<:Real}
    error("Chosen gradient order not supported.")
end