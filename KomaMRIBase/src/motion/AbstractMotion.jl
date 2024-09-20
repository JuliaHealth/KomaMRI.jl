abstract type AbstractMotion{T<:Real} end

# NoMotion
include("nomotion/NoMotion.jl")

# MotionList
include("motionlist/Action.jl")
include("motionlist/SpinSpan.jl")
include("motionlist/TimeSpan.jl")
include("motionlist/Motion.jl")
include("motionlist/MotionList.jl")