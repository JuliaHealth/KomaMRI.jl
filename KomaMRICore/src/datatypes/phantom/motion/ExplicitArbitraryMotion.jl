@with_kw mutable struct ExplicitArbitraryMotion{T} <: MotionModel{T}
	Ux   ::Union{Nothing, AbstractArray{T, 2}} = nothing
    Uy   ::Union{Nothing, AbstractArray{T, 2}} = nothing
    Uz   ::Union{Nothing, AbstractArray{T, 2}} = nothing
    flags::Union{Nothing, AbstractArray{T, 2}} = nothing #Boolean
end


Base.getindex(motion::ExplicitArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                               q::Union{AbstractRange,AbstractVector,Colon}) = begin
    ExplicitArbitraryMotion(Ux    = motion.Ux    !== nothing ? motion.Ux[p,q]    : nothing,
                            Uy    = motion.Uy    !== nothing ? motion.Uy[p,q]    : nothing,
                            Uz    = motion.Uz    !== nothing ? motion.Uz[p,q]    : nothing,
                            flags = motion.flags !== nothing ? motion.flags[p,q] : nothing)
end


function get_displacements(motion::ExplicitArbitraryMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}
    ex = (length(t)==1)
    
    Ux = motion.Ux    !== nothing ? (ex ? reshape(motion.Ux,(length(x),)) : hcat(motion.Ux, motion.Ux[:,end]))  : nothing
    Uy = motion.Uy    !== nothing ? (ex ? reshape(motion.Uy,(length(x),)) : hcat(motion.Uy, motion.Uy[:,end]))  : nothing
    Uz = motion.Uz    !== nothing ? (ex ? reshape(motion.Uz,(length(x),)) : hcat(motion.Uz, motion.Uz[:,end]))  : nothing
    Ux, Uy, Uz, motion.flags
end
