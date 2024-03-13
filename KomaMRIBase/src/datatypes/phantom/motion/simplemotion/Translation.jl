"""
Translation Motion

Parameters:
- Initial time      (ti)
- Final time        (tf)
- Translation in x  (dx)
- Translation in y  (dy)
- Translation in z  (dz)

     [ 0,                   t <= ti
ux = [ dx*(t-ti)/(tf-ti),   ti < t < tf
     [ dx,                  t >= tf

     [ 0,                   t <= ti
uy = [ dy*(t-ti)/(tf-ti),   ti < t < tf
     [ dy,                  t >= tf

     [ 0,                   t <= ti
uz = [ dz*(t-ti)/(tf-ti),   ti < t < tf
     [ dz,                  t >= tf


"""

@with_kw struct Translation{T<:Real} <: SimpleMotionType{T}
    ti::T = 0.0
    tf::T = 0.0
    dx::T = 0.0
    dy::T = 0.0
    dz::T = 0.0
end

displacement_x(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    (t .<= motion_type.ti) .* 0    .+    ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.dx .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+    (t .>= motion_type.tf) .* motion_type.dx
end

displacement_y(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    (t .<= motion_type.ti) .* 0    .+    ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.dy .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+    (t .>= motion_type.tf) .* motion_type.dy
end

displacement_z(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    (t .<= motion_type.ti) .* 0    .+    ((t .> motion_type.ti) .& (t .< motion_type.tf)) .* motion_type.dz .* (t .- motion_type.ti) ./ (motion_type.tf - motion_type.ti)   .+    (t .>= motion_type.tf) .* motion_type.dz
end
