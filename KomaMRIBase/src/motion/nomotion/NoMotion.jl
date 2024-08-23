struct NoMotion{T<:Real} <: AbstractMotionSet{T} end

Base.getindex(mv::NoMotion, p::AbstractVector) = mv 
Base.view(mv::NoMotion, p::AbstractVector)     = mv

""" Addition of NoMotions """
Base.vcat(m1::NoMotion{T}, m2::NoMotion{T}, Ns1::Int, Ns2::Int) where {T<:Real} = NoMotion{T}()
function Base.vcat(m1::NoMotion{T}, m2::AbstractMotionSet{T}, Ns1::Int, Ns2::Int) where {T<:Real}
    mv_aux = Motion{T}[]
    for m in m2.motions
        m_aux = copy(m)
        if m_aux.spins == Colon()
            m_aux.spins = 1:Ns2
        end
        m_aux.spins = m_aux.spins .+ Ns1
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux)
end
function Base.vcat(m1::AbstractMotionSet{T}, m2::NoMotion{T}, Ns1::Int, Ns2::Int) where {T<:Real}
    mv_aux = Motion{T}[]
    for m in m1.motions
        m_aux = copy(m)
        if m_aux.spins == Colon()
            m_aux.spins = 1:Ns1
        end
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux)
end

Base.:(==)(m1::NoMotion{T}, m2::NoMotion{T}) where {T<:Real} = true
Base.:(≈)(m1::NoMotion{T}, m2::NoMotion{T}) where {T<:Real}  = true

function get_spin_coords(
    mv::NoMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    return x, y, z
end

"""
    times = times(motion)
"""
times(mv::NoMotion{T}) where {T<:Real} = [zero(T)]

"""
    sort_motions!
"""
sort_motions!(mv::NoMotion) =  nothing