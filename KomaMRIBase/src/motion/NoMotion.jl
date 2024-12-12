"""
    nomotion = NoMotion()

NoMotion struct. It is used to create static phantoms.

# Returns
- `nomotion`: (`::NoMotion`) NoMotion struct

# Examples
```julia-repl
julia> nomotion = NoMotion()
```
"""
struct NoMotion end

Base.getindex(mv::NoMotion, p) = mv 
Base.view(mv::NoMotion, p)     = mv

""" Addition of NoMotions """
# NoMotion + NoMotion
Base.vcat(m1::NoMotion,   m2::NoMotion, Ns1, Ns2) = m1
# NoMotion + MotionList
Base.vcat(m1::MotionList, m2::NoMotion, Ns1, Ns2) = vcat(m2, m1, 0, Ns1)
function Base.vcat(m1::NoMotion, m2::MotionList{T}, Ns1, Ns2) where {T}
    mv_aux = Motion{T}[]
    for m in m2.motions
        m_aux = copy(m)
        m_aux.spins = expand(m_aux.spins, Ns2)
        m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux)
end
# NoMotion + Motion
Base.vcat(m1::Motion, m2::NoMotion, Ns1, Ns2) = vcat(m2, m1, 0, Ns1)
function Base.vcat(m1::NoMotion, m2::Motion{T}, Ns1, Ns2) where {T}
    m_aux = copy(m2)
    m_aux.spins = expand(m_aux.spins, Ns2)
    m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
    return m_aux
end

""" Compare two NoMotions """
Base.:(==)(m1::NoMotion, m2::NoMotion) = true
Base.:(≈)(m1::NoMotion, m2::NoMotion)   = true

function get_spin_coords(
    mv::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t
) where {T<:Real}
    return x, y, z
end
add_jump_times!(t, ::NoMotion) = nothing