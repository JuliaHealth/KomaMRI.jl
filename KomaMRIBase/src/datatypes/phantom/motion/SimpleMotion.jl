@with_kw mutable struct SimpleMotion{T} <: MotionModel{T}
    ux::Function = (x,y,z,t)->0
	uy::Function = (x,y,z,t)->0
	uz::Function = (x,y,z,t)->0
end

Base.getindex(motion::SimpleMotion, p::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractRange, q::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractVector) = motion

function get_U(motion::SimpleMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}                      
    Ns = length(x)
                        
    Ux = zeros(Ns,1) .+ motion.ux(x, y, z, t')
    Uy = zeros(Ns,1) .+ motion.uy(x, y, z, t')
    Uz = zeros(Ns,1) .+ motion.uz(x, y, z, t')

    Ux = sum(abs.(Ux);dims=2) != zeros(Ns,1) ? Ux : nothing
    Uy = sum(abs.(Uy);dims=2) != zeros(Ns,1) ? Uy : nothing
    Uz = sum(abs.(Uz);dims=2) != zeros(Ns,1) ? Uz : nothing

    Ux,Uy,Uz
end


"""
    motion = initialize_motion(obj.motion, seqd.t)
"""
function initialize_motion(motion::SimpleMotion, t::AbstractVector{T}, sim_params::Dict) where {T<:Real}
    return motion
end


function is_dynamic(motion::SimpleMotion)
    x = 0:0.01:10
    y = 0:0.01:10
    z = 0:0.01:10
    t = 0:0.01:10

    Ux,Uy,Uz = get_U(motion,x,y,z,t)
    return reduce(|,([Ux,Uy,Uz] .!== nothing))
end


function get_displacements(motion::SimpleMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}
	Ux = motion.ux(x, y, z, t')
    Uy = motion.uy(x, y, z, t')
    Uz = motion.uz(x, y, z, t')

    Ux, Uy, Uz, nothing
end