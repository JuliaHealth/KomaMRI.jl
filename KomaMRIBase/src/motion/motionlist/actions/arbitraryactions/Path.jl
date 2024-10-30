@doc raw"""
    path = Path(dx, dy, dz)

Path struct. For this action (and for `FlowPath`),
motion is not defined solely on the basis of 
three numerical parameters, one for each spatial direction,
as occurs for the `Translate`, `Rotate` and `HeartBeat` actions.

For this action, it is necessary to define 
motion for each spin independently, in x (`dx`), y (`dy`) and z (`dz`).
`dx`, `dy` and `dz` are now three matrixes, of (``N_{spins}* \times \; N_{discrete\,times}``) each.
This means that each row corresponds to a spin trajectory over a set of discrete time instants.

!!! note
    *When creating a motion with `Flow` or `FlowPath`, you must make sure that 
    the number of rows of the matrices `dx`, `dy` and `dz` matches the number 
    of spins that are affected by the motion. 

    Remember that the range of spins affected by a motion 
    is defined by the `spins` field of the `Motion` struct
    
    example:
    ```julia-repl
    julia> motion = Motion(
        action = Path(
            dx=[0.01 0.02; 0.02 0.03],  # 2 rows
            dy=[0.02 0.03; 0.03 0.04], 
            dz=[0.03 0.04; 0.04 0.05]),
        time = TimeRange(0.0, 1.0),
        spins = SpinRange(1:2)          # 2 spins
    )
    ```

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z

# Returns
- `path`: (`::Path`) Path struct

# Examples
```julia-repl
julia> path = Path(
           dx=[0.01 0.02; 0.02 0.03], 
           dy=[0.02 0.03; 0.03 0.04], 
           dz=[0.03 0.04; 0.04 0.05]
       )
```
"""
@with_kw struct Path{T<:Real} <: ArbitraryAction{T}
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
end