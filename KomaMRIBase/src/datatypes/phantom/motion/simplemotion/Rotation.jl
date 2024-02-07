"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

@with_kw mutable struct Rotation{T<:Real} <: SimpleMotionType
    axis::AbstractVector{T}     # Rotation axis vector
    point::AbstractVector{T}    # Rotation axis point
    f::T                        # Angular velocity [Hz]
end

function SimpleMotion(rot::Rotation)
    # Rotation matrix
    R(t) = axis_angle.(Ref(rot.axis), 2π*rot.f*t)

    # Displace using the rotation point
    displace(x,y,z) =  [[x,y,z][i] .- rot.point[i] for i in 1:length(rot.point)]

    # Conversion of x,y,z into a (3 x Ns) matrix
    matrix(x,y,z) = hcat(x,y,z)'

    rotation(x,y,z,t) = begin
        xd, yd, zd = displace(x,y,z)
        xyz = matrix(xd,yd,zd)
        display(R(t))
        return cat(map(x -> x*xyz, R(t))...,dims=3)
    end

    ux(x,y,z,t) = @view(rotation(x,y,z,t)[1,:,:]) .+ (rot.point[1] .- x)  
    uy(x,y,z,t) = @view(rotation(x,y,z,t)[2,:,:]) .+ (rot.point[2] .- y)  
    uz(x,y,z,t) = @view(rotation(x,y,z,t)[3,:,:]) .+ (rot.point[3] .- z)  

    return SimpleMotion(type=rot,ux=ux,uy=uy,uz=uz)
end