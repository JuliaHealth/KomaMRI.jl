struct FlowTrajectory{T<:Real, TS<:TimeScale{T}} <: ArbitraryMotion{T}
    times::TS
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    resetmag::AbstractArray{Bool}
end