struct NoMotion <: MotionModel end

Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                q::Union{AbstractRange,AbstractVector,Colon}) = motion

displacement_x(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0
displacement_y(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0
displacement_z(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0