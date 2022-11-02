"""
    y = trapz(Δt, x)

Trapezoidal integration for every spin of a phantom.

!!! note
    In practice, this function is used to integrate (Gx * x + Gy * y + Gz * z) * Δt for all
    the spins. `NΔt` is the length of `Δt`. `Ns` stands for the number of spins of a phantom.
    `x` is a matrix which rows represents different spins and columns are different times
    and the elements are the field Gx * x + Gy * y + Gz * z values.

# Arguments
- `Δt`: (`1 x NΔt ::Matrix{Float64}`, `[s]`) delta time 1-row array
- `x`: (`Ns x (NΔt+1) ::Matrix{Float64}`, `[T]`) magnitude of the field Gx * x + Gy * y +
    Gz * z

# Returns
- `y`: (`Ns x 1 ::Matrix{Float64}`, `[T*s]`) vector where every element is the integral
    of (Gx * x + Gy * y + Gz * z) * Δt for every spin of a phantom
"""
function trapz(Δt, x)
    y = 1/2 * (x[:, 2:end] .+ x[:, 1:end-1]) .* Δt
    y = sum(y, dims=2)
    return y
end

"""
    y = cumtrapz(Δt, x)

Trapezoidal cumulative integration over time for every spin of a phantom.

!!! note
    In practice, this function is used to get the cumulative integration of
    (Gx * x + Gy * y + Gz * z) * Δt for all the spins. `NΔt` is the length of `Δt`. `Ns`
    stands for the number of spins of a phantom. `x` is a matrix which rows represents
    different spins and columns are different times and the elements are the field Gx * x +
    Gy * y + Gz * z values. Note that the last column of the result of this function
    [`cumtrapz`](@ref) is equal to the result of the function [`trapz`](@ref)

# Arguments
- `Δt`: (`1 x NΔt ::Matrix{Float64}`, `[s]`) delta time 1-row array
- `x`: (`Ns x (NΔt+1) ::Matrix{Float64}`, `[T]`) magnitude of the field Gx * x + Gy * y +
    Gz * z

# Returns
- `y`: (`Ns x NΔt ::Matrix{Float64}`, `[T*s]`) matrix where every column is the
    cumulative integration over time of (Gx * x + Gy * y + Gz * z) * Δt for every spin of a
    phantom
"""
function cumtrapz(Δt, x)
    y = 1/2 * (x[:, 2:end] .+ x[:, 1:end-1]) .* Δt
    y = cumsum(y, dims=2)
    return y
end
