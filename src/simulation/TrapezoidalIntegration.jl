"""
Trapezoidal integration of x Nx1 and Δt 1xM to obtain y NxM.
"""
function trapz(Δt,x)
    y = 1/2 * (x[:,2:end] .+ x[:,1:end-1]) .* Δt
    y = sum(y, dims=2)
end
"""
Trapezoidal cumulative integration of x Nx1 and Δt 1xM to obtain y NxM.
"""
function cumtrapz(Δt,x)
    y = 1/2 * (x[:,2:end] .+ x[:,1:end-1]) .* Δt
    y = cumsum(y, dims=2)
end