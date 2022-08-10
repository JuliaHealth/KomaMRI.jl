"""
    trapz(Δt, x)

Trapezoidal integration.

# Arguments
- `Δt`: dimension 1xM
- `x`: dimension Nx1

# Returns
- `y`: dimension NxM
"""
function trapz(Δt, x)
    y = 1/2 * (x[:, 2:end] .+ x[:, 1:end-1]) .* Δt
    y = sum(y, dims=2)
    return y
end

"""
    cumtrapz(Δt, x)

Trapezoidal cumulative integration.

# Arguments
- `Δt`: dimension 1xM
- `x`: dimension Nx1

# Returns
- `y`: dimension NxM
"""
function cumtrapz(Δt, x)
    y = 1/2 * (x[:, 2:end] .+ x[:, 1:end-1]) .* Δt
    y = cumsum(y, dims=2)
    return y
end
