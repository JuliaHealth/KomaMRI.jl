struct FlowTrajectory{T<:Real, TS<:AbstractTimeSpan{T}} <: ArbitraryMotion{T}
    time::TS
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    spin_reset::AbstractArray{Bool}
end