abstract type AbstractMotion{T<:Real} end

is_composable(m::AbstractMotion) = false

struct MotionList{T<:Real} <: AbstractMotionList{T}
    motions::Vector{<:AbstractMotion{T}}
end

MotionList(motions...) = length([motions]) > 0 ? MotionList([motions...]) : @error "You must provide at least one motion as input argument. If you do not want to define motion, use `NoMotion{T}()`"

include("motion/SimpleMotion.jl")
include("motion/ArbitraryMotion.jl")

Base.getindex(mv::MotionList, p::Union{AbstractRange, AbstractVector, Colon, Integer}) = MotionList(getindex.(mv.motions, Ref(p))) 
Base.view(mv::MotionList, p::Union{AbstractRange, AbstractVector, Colon, Integer})     = MotionList(view.(mv.motions, Ref(p)))

""" Addition of MotionLists """
function Base.vcat(m1::MotionList{T}, m2::MotionList{T}, Ns1::Int, Ns2::Int) where {T<:Real}
    mv1 = m1.motions
    mv1_aux = AbstractMotion{T}[]
    for i in 1:length(mv1)
        if typeof(mv1[i]) <: ArbitraryMotion
            zeros1 = similar(mv1[i].dx, Ns2, size(mv1[i].dx, 2))
            zeros1 .= zero(T)
            push!(mv1_aux, typeof(mv1[i])(mv1[i].time, [[getfield(mv1[i], d); zeros1] for d in filter(x -> x != :time, fieldnames(typeof(mv1[i])))]...)) 
        else
            push!(mv1_aux, mv1[i])
        end
    end
    mv2 = m2.motions
    mv2_aux = AbstractMotion{T}[]
    for i in 1:length(mv2)
        if typeof(mv2[i]) <: ArbitraryMotion
            zeros2 = similar(mv2[i].dx, Ns1, size(mv2[i].dx, 2))
            zeros2 .= zero(T)
            push!(mv2_aux, typeof(mv2[i])(mv2[i].time, [[zeros2; getfield(mv2[i], d)] for d in filter(x -> x != :time, fieldnames(typeof(mv2[i])))]...))
        else
            push!(mv2_aux, mv2[i])
        end
    end
    return MotionList([mv1_aux; mv2_aux])
end

""" Compare two motion vectors """
function Base.:(==)(mv1::MotionList{T}, mv2::MotionList{T}) where {T<:Real}
    sort_motions!(mv1)
    sort_motions!(mv2)
    return reduce(&, mv1.motions .== mv2.motions)
end
function Base.:(≈)(mv1::MotionList{T}, mv2::MotionList{T}) where {T<:Real} 
    sort_motions!(mv1)
    sort_motions!(mv2)
    return reduce(&, mv1.motions .≈ mv2.motions)
end

"""
    x, y, z = get_spin_coords(motion, x, y, z, t)

Calculates the position of each spin at a set of arbitrary time instants, i.e. the time steps of the simulation. 
For each dimension (x, y, z), the output matrix has ``N_{\text{spins}}`` rows and `length(t)` columns.

# Arguments
- `motion`: (`::Vector{<:AbstractMotion{T<:Real}}`) phantom motion
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `t`: (`::AbstractArray{T<:Real}`) horizontal array of time instants

# Returns
- `x, y, z`: (`::Tuple{AbstractArray, AbstractArray, AbstractArray}`) spin positions over time
"""
function get_spin_coords(
    mv::MotionList{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    # Buffers for positions:
    xt, yt, zt = x .+ 0*t, y .+ 0*t, z .+ 0*t
    # Buffers for displacements:
    ux, uy, uz = similar(xt), similar(yt), similar(zt) 

    # Composable motions: they need to be run sequentially. Note that they depend on xt, yt , and zt
    for m in Iterators.filter(is_composable, mv.motions)
        displacement_x!(ux, m, xt, yt, zt, t)
        displacement_y!(uy, m, xt, yt, zt, t)
        displacement_z!(uz, m, xt, yt, zt, t)
        xt .+= ux
        yt .+= uy
        zt .+= uz
    end
    # Additive motions: these motions can be run in parallel
    for m in Iterators.filter(!is_composable, mv.motions)
        displacement_x!(ux, m, x, y, z, t)
        displacement_y!(uy, m, x, y, z, t)
        displacement_z!(uz, m, x, y, z, t)
        xt .+= ux
        yt .+= uy
        zt .+= uz
    end
    return xt, yt, zt
end

"""
    times = times(motion)
"""
times(m::AbstractMotion) = times(m.time)
times(mv::MotionList{T}) where {T<:Real} = begin
    nodes = reduce(vcat, [times(m) for m in mv.motions]; init=[zero(T)])
    return unique(sort(nodes))
end

"""
    sort_motions!(motion_list)
sort_motions motions in a list according to their starting time
"""
function sort_motions!(mv::MotionList{T}) where {T<:Real}
    sort!(mv.motions; by=m -> times(m)[1])
    return nothing
end