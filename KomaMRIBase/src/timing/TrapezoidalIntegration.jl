"""
    y = trapz(Δt, x)

Trapezoidal integration. In practice, this function integrates the expression
(Gx * x + Gy * y + Gz * z) * Δt over time for all spins. `NΔt` represents the length of the
time interval `Δt`, and `Ns` denotes the number of spins in a phantom. The matrix `x` has
rows representing different spins and columns representing different times, with elements
corresponding to the field values (Gx * x + Gy * y + Gz * z).

# Arguments
- `Δt`: (`1 x NΔt ::AbstractArray{Real}`, `[s]`) delta time 1-row array
- `x`: (`Ns x (NΔt+1) ::AbstractArray{Real}`, `[T]`) magnitude of the field
    Gx * x + Gy * y + Gz * z

# Returns
- `y`: (`Ns x 1 ::AbstractArray{Real}`, `[T*s]`) vector where every element is the integral
    of (Gx * x + Gy * y + Gz * z) * Δt for every spin of a phantom
"""
function trapz(Δt::AbstractArray{T}, x::AbstractArray{T}) where {T<:Real}
    y = (x[:, 2:end] .+ x[:, 1:end-1]) .* (Δt / 2)
    y = sum(y, dims=2)
    return y
end

"""
    y = cumtrapz(Δt, x)

Trapezoidal cumulative integration. Same as [`KomaMRIBase.trapz`](@ref) but it store
cumulative integration results.

# Arguments
- `Δt`: (`1 x NΔt ::AbstractArray{Real}`, `[s]`) delta time 1-row array
- `x`: (`Ns x (NΔt+1) ::AbstractArray{Real}`, `[T]`) magnitude of the field Gx * x + Gy * y +
    Gz * z

# Returns
- `y`: (`Ns x NΔt ::AbstractArray{Real}`, `[T*s]`) matrix where every column is the
    cumulative integration over time of (Gx * x + Gy * y + Gz * z) * Δt for every spin of a
    phantom
"""
function cumtrapz(Δt::AbstractArray{T}, x::AbstractArray{T}) where {T<:Real}
    y =  (x[:, 2:end] .+ x[:, 1:end-1]) .* (Δt / 2)
    y = cumsum(y, dims=2)
    return y
end
