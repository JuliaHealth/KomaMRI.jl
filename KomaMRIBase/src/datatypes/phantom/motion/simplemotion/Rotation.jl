"""
Rotation Motion
 
Parameters:
- Initial time   (ti)
- Final time     (tf)
- Rotation in x  (pitch)
- Rotation in y  (roll)
- Rotation in z  (yaw)

        [ 0,                      t <= ti
pitch = [ pitch*(t-ti)/(tf-ti),   ti < t < tf
        [ pitch,                  t >= tf

        [ 0,                      t <= ti
roll =  [ roll*(t-ti)/(tf-ti),    ti < t < tf
        [ roll,                   t >= tf

        [ 0,                      t <= ti
yaw =   [ yaw*(t-ti)/(tf-ti),     ti < t < tf
        [ yaw,                    t >= tf


x' = cos(yaw)cos(roll)x + (cos(yaw)sin(roll)sin(pitch) - sin(yaw)cos(pitch))y + (cos(yaw)sin(roll)cos(pitch) + sin(yaw)sin(pitch))z
y' = sin(yaw)cos(roll)x + (sin(yaw)sin(roll)sin(pitch) + cos(yaw)cos(pitch))y + (sin(yaw)sin(roll)cos(pitch) - cos(yaw)sin(pitch))z
z' = -sin(roll)x        + cos(roll)sin(pitch)y                                + cos(roll)cos(pitch)z

ux = x' - x
uy = y' - y
uz = z' - z

"""

@with_kw struct Rotation{T<:Real} <: SimpleMotionType{T}
    ti::T = 0.0
    tf::T = 0.0
    pitch::T = 0.0
    roll::T = 0.0
    yaw::T = 0.0
end 

displacement_x(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    pitch   = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.pitch  .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.pitch
    roll    = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.roll   .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.roll
    yaw     = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.yaw    .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.yaw  
    return cos.(yaw) .* cos.(roll) .* x   +   (cos.(yaw) .* sin.(roll) .* sin.(pitch) .- sin.(yaw) .* cos.(pitch)) .* y   +   (cos.(yaw) .* sin.(roll) .* cos.(pitch) .+ sin.(yaw) .* sin.(pitch)) .*z   .-   x
end

displacement_y(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    pitch   = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.pitch  .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.pitch
    roll    = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.roll   .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.roll
    yaw     = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.yaw    .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.yaw  
    return sin.(yaw) .* cos.(roll) .* x   +   (sin.(yaw) .* sin.(roll) .* sin.(pitch) .+ cos.(yaw) .* cos.(pitch)) .* y   +   (sin.(yaw) .* sin.(roll) .* cos.(pitch) .- cos.(yaw) .* sin.(pitch)) .* z  .-  y
end

displacement_z(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    pitch   = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.pitch  .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.pitch
    roll    = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.roll   .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.roll
    yaw     = (t .<= motion_type.ti) .* 0   .+   ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.yaw    .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+   (t .>= motion_type.tf) .* motion_type.yaw  
    return -sin.(roll) .* x   +   cos.(roll) .* sin.(pitch) .* y  .-  z 
end