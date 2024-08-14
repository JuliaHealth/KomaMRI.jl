@doc raw"""
    flowtrajectory = FlowTrajectory(time, dx, dy, dz)

FlowTrajectory motion struct. (...)

# Arguments
- `time`: (`::AbstractTimeSpan{T<:Real}`, `[s]`) time scale
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags

# Returns
- `flowtrajectory`: (`::FlowTrajectory`) FlowTrajectory struct

# Examples
```julia-repl
julia> ftr = FlowTrajectory(time=TimeRange(0.0, 0.5), dx=[0.01 0.02], dy=[0.02 0.03], dz=[0.03 0.04], spin_reset=[false, false])
```
"""
struct FlowTrajectory{T<:Real, TS<:AbstractTimeSpan{T}} <: ArbitraryMotion{T}
    time::TS
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    spin_reset::AbstractArray{Bool}
end