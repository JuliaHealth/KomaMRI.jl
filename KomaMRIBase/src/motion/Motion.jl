"""
    motion = Motion(action)
    motion = Motion(action, time)
    motion = Motion(action, time, spins)

Motion struct. It defines the motion, during a certain time interval,
of a given group of spins. It is composed by three fields: `action`, which 
defines the motion itself, `time`, which accounts for the time during
which the motion takes place, and `spins`, which indicates the spins 
that are affected by that motion.

# Arguments
- `action`: (`::AbstractAction{T<:Real}`) action, such as [`Translate`](@ref) or [`Rotate`](@ref)
- `time`: (`::TimeCurve{T<:Real}`, `=TimeRange(0.0)`) time information about the motion
- `spins`: (`::AbstractSpinSpan`, `=AllSpins()`) spin indexes affected by the motion

# Returns
- `motion`: (`::Motion`) Motion struct

# Examples
```julia-repl
julia> motion =  Motion(
            action = Translate(0.01, 0.0, 0.02),
            time = TimeRange(0.0, 1.0),
            spins = SpinRange(1:10)
       )
```
"""
@with_kw mutable struct Motion{T<:Real}
    action::AbstractAction{T}
    time  ::TimeCurve{T}      = TimeRange(t_start=zero(typeof(action).parameters[1]), t_end=eps(typeof(action).parameters[1]))
    spins ::AbstractSpinSpan  = AllSpins()
end

# Main constructors
function Motion(action) 
    T = first(typeof(action).parameters)
    return Motion(action, TimeRange(t_start=zero(T), t_end=eps(T)), AllSpins())
end
function Motion(action, time::TimeCurve)
    T = first(typeof(action).parameters)
    return Motion(action, time, AllSpins())
end
function Motion(action, spins::AbstractSpinSpan)
    T = first(typeof(action).parameters)
    return Motion(action, TimeRange(t_start=zero(T), t_end=eps(T)), spins)
end

# Custom constructors
"""
    tr = translate(dx, dy, dz, time, spins)

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `tr`: (`::Motion`) Motion struct

# Examples
```julia-repl
julia> tr = translate(0.01, 0.02, 0.03, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function translate(dx, dy, dz, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(Translate(dx, dy, dz), time, spins)
end

"""
    rt = rotate(pitch, roll, yaw, spins)

# Arguments
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Keywords
- `center`: (`::NTuple{3,Real}` or `::CenterOfMass`) center of rotation, given in global coordinates. Default is center of mass.

# Returns
- `rt`: (`::Motion`) Motion struct with [`Rotate`](@ref) action

# Examples
```julia-repl
julia> rt = rotate(15.0, 0.0, 20.0, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function rotate(pitch, roll, yaw, time=TimeRange(t_start=zero(eltype(pitch)), t_end=eps(eltype(pitch))), spins=AllSpins(); center=CenterOfMass())
    return Motion(Rotate(pitch, roll, yaw, center), time, spins)
end

"""
    hb = heartbeat(circumferential_strain, radial_strain, longitudinal_strainl, time, spins)

# Arguments
- `circumferential_strain`: (`::Real`) contraction parameter
- `radial_strain`: (`::Real`) contraction parameter
- `longitudinal_strain`: (`::Real`) contraction parameter
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `hb`: (`::Motion`) Motion struct with [`HeartBeat`](@ref) action

# Examples
```julia-repl
julia> hb = heartbeat(-0.3, -0.2, 0.0, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function heartbeat(circumferential_strain, radial_strain, longitudinal_strain, time=TimeRange(t_start=zero(eltype(circumferential_strain)), t_end=eps(eltype(circumferential_strain))), spins=AllSpins())
    return Motion(HeartBeat(circumferential_strain, radial_strain, longitudinal_strain), time, spins)
end

"""
    pt = path(dx, dy, dz, time, spins)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `pt`: (`::Motion`) Motion struct with [`Path`](@ref) action

# Examples
```julia-repl
julia> pt = path(
          [0.01 0.02; 0.02 0.03], 
          [0.02 0.03; 0.03 0.04], 
          [0.03 0.04; 0.04 0.05], 
          TimeRange(0.0, 1.0), 
          SpinRange(1:10)
       )
```
"""
function path(dx, dy, dz, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(Path(dx, dy, dz), time, spins)
end

"""
    fp = flowpath(dx, dy, dz, spin_reset, time, spins)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `fp`: (`::Motion`) Motion struct with [`FlowPath`](@ref) action

# Examples
```julia-repl
julia> fp = flowpath(
          [0.01 0.02; 0.02 0.03], 
          [0.02 0.03; 0.03 0.04], 
          [0.03 0.04; 0.04 0.05], 
          [false false; false true],
          TimeRange(0.0, 1.0), 
          SpinRange(1:10)
       )
```
"""
function flowpath(dx, dy, dz, spin_reset, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(FlowPath(dx, dy, dz, spin_reset), time, spins)
end

""" Compare two Motions """
Base.:(==)(m1::Motion, m2::Motion) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::Motion, m2::Motion)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

""" Motion sub-group """
function Base.getindex(m::Motion, p)
    idx, spin_range = m.spins[p]
    return idx !== nothing ? Motion(m.action[idx], m.time, spin_range) : NoMotion()
end
function Base.view(m::Motion, p)
    idx, spin_range = @view(m.spins[p])
    return idx !== nothing ? Motion(@view(m.action[idx]), m.time, spin_range) : NoMotion()
end

"""
    x, y, z = get_spin_coords(motion, x, y, z, t)

Calculates the position of each spin at a set of arbitrary time instants, i.e. the time steps of the simulation. 
For each dimension (x, y, z), the output matrix has ``N_{\t{spins}}`` rows and `length(t)` columns.

# Arguments
- `motion`: (`::Union{NoMotion, Motion{T<:Real} MotionList{T<:Real}}`) phantom motion
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `t`: horizontal array of time instants

# Returns
- `x, y, z`: (`::Tuple{AbstractArray, AbstractArray, AbstractArray}`) spin positions over time
"""
function get_spin_coords(
    m::Motion{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t
) where {T<:Real}
    ux, uy, uz = x .* (0*t), y .* (0*t), z .* (0*t) # Buffers for displacements
    t_unit = unit_time(t, m.time)
    idx = get_indexing_range(m.spins)
    displacement_x!(@view(ux[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    displacement_y!(@view(uy[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    displacement_z!(@view(uz[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    return x .+ ux, y .+ uy, z .+ uz
end

# Auxiliary functions
times(m::Motion) = times(m.time)
is_composable(m::Motion) = is_composable(m.action)

"""
    add_key_time_points!(t, motion)
"""
function add_key_time_points!(t, m::Motion)
    add_key_time_points!(t, m.action, m.time.t_start, m.time.t_end, m.time.periods, m.time.periodic)
    return nothing
end
function add_key_time_points!(t, a, t_start::T, t_end::T, periods, periodic) where T
    isempty(t) && return
    aux = T[] 
    period = sum((t_end - t_start) .* periods)
    t_max = maximum(t)
    add_period_times!(aux, t_start, t_end, periods)
    add_reset_times!(aux, a, t_start, t_end, periods)
    extend_periodic!(aux, t_max, period, Val(periodic))
    append!(t, aux[aux .<= t_max])
    return nothing
end

"""
    extend_periodic!(aux, t_max, period, periodic)
"""
function extend_periodic!(aux, t_max, period, periodic::Val{false})
    return nothing
end
function extend_periodic!(aux, t_max, period, periodic::Val{true})
    n_periods = floor(Int, t_max / period)
    if n_periods > 0
        initial_size = length(aux)
        sizehint!(aux, initial_size * (n_periods + 1))
        for n in 1:n_periods
            append!(aux, aux[1:initial_size] .+ n*period)
        end
    end
    return nothing
end

"""
    add_period_times!(t, t_start, t_end, periods)
"""
function add_period_times!(t, t_start, t_end, periods)
    period_times = times([t_start, t_end], t_start, t_end, periods)
    append!(t, period_times .+ MIN_RISE_TIME .* ((-1) .^ ((1:length(period_times)) .+ 1)))
    return nothing
end

"""
    add_reset_times!(t, action, t_start, t_end, periods)
""" 
function add_reset_times!(t, ::AbstractAction, t_start, t_end, periods)
    return nothing 
end