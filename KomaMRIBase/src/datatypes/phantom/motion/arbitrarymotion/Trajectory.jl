@doc raw"""
    trajectory = Trajectory(time, dx, dy, dz)

Trajectory motion struct. (...)

# Arguments
- `time`: (`::AbstractTimeSpan{T<:Real}`, `[s]`) time scale
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z

# Returns
- `trajectory`: (`::Trajectory`) Trajectory struct

# Examples
```julia-repl
julia> tr = Trajectory(time=TimeRange(0.0, 0.5), dx=[0.01 0.02], dy=[0.02 0.03], dz=[0.03 0.04])
```
"""
struct Trajectory{T<:Real, TS<:AbstractTimeSpan{T}} <: ArbitraryMotion{T}
    time::TS
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
end