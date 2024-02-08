"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

mutable struct Rotation{T<:Real} <: SimpleMotionType
    axis::AbstractVector{T}     # Rotation axis vector
    point::AbstractVector{T}    # Rotation axis point
    f::T                        # Angular velocity [Hz]
end

function SimpleMotion(type::Rotation)
    # Rotation matrix
    # R(t) = axis_angle.(Ref(type.axis), 2π*type.f*t)
    R(t) = [axis_angle(type.axis,x) for x in 2π*type.f*t] # Scalar indexing. Solve

    # Displace using the rotation point
    displace(x,y,z) =  [[x,y,z][i] .- type.point[i] for i in 1:length(type.point)]

    # Conversion of x,y,z into a (3 x Ns) matrix
    matrix(x,y,z) = hcat(x,y,z)'

    rotation(x,y,z,t) = begin
        xd, yd, zd = displace(x,y,z)
        xyz = matrix(xd,yd,zd)
        display(R(t))
        rotated = map(x -> x*xyz, R(t))
        return cat(rotated...,dims=3)
    end

    ux(x,y,z,t) = @view(rotation(x,y,z,t)[1,:,:]) .+ (type.point[1] .- x)  
    uy(x,y,z,t) = @view(rotation(x,y,z,t)[2,:,:]) .+ (type.point[2] .- y)  
    uz(x,y,z,t) = @view(rotation(x,y,z,t)[3,:,:]) .+ (type.point[3] .- z)  

    return SimpleMotion(type,ux,uy,uz)
end