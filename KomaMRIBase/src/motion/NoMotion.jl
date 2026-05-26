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
function Base.vcat(m1::NoMotion, m2::MotionList, Ns1, Ns2)
    mv_aux = Motion[]
    for m in m2.motions
        m_aux = deepcopy(m)
        m_aux.spins = expand(m_aux.spins, Ns2)
        m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux)
end
# NoMotion + Motion
Base.vcat(m1::Motion, m2::NoMotion, Ns1, Ns2) = vcat(m2, m1, 0, Ns1)
function Base.vcat(m1::NoMotion, m2::Motion, Ns1, Ns2)
    m_aux = deepcopy(m2)
    m_aux.spins = expand(m_aux.spins, Ns2)
    m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
    return m_aux
end

""" Compare two NoMotions """
Base.:(==)(m1::NoMotion, m2::NoMotion) = true
Base.:(≈)(m1::NoMotion, m2::NoMotion)   = true

get_spin_coords(::NoMotion, x::AbstractVector, y::AbstractVector, z::AbstractVector, t) = (x, y, z)
add_key_time_points!(t, ::NoMotion) = nothing
