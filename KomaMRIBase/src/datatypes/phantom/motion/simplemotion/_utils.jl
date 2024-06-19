"""
    sorted_motion = sort_motion(motion)

Sorts, with respect to time, the motion types of a `SimpleMotion` instance.
No allocations, since it uses the TupleTools.sort method
"""
function sort_motion(motion::SimpleMotion)
    return SimpleMotion(sort(motion.types; by=m -> times(m)[1])) 
end


"""
    sort(t::Tuple; lt=isless, by=identity, rev::Bool=false) -> ::Tuple

Sorts the tuple `t`. Extracted from TupleTools.jl
"""
sort(t::Tuple; lt=isless, by=identity, rev::Bool=false) = _sort(t, lt, by, rev)
@inline function _sort(t::Tuple, lt=isless, by=identity, rev::Bool=false)
    t1, t2 = _split(t)
    t1s = _sort(t1, lt, by, rev)
    t2s = _sort(t2, lt, by, rev)
    return _merge(t1s, t2s, lt, by, rev)
end
_sort(t::Tuple{Any}, lt=isless, by=identity, rev::Bool=false) = t
_sort(t::Tuple{}, lt=isless, by=identity, rev::Bool=false) = t

function _split(t::Tuple)
    N = length(t)
    M = N >> 1
    return ntuple(i -> t[i], M), ntuple(i -> t[i + M], N - M)
end

function _merge(t1::Tuple, t2::Tuple, lt, by, rev)
    if lt(by(first(t1)), by(first(t2))) != rev
        return (first(t1), _merge(tail(t1), t2, lt, by, rev)...)
    else
        return (first(t2), _merge(t1, tail(t2), lt, by, rev)...)
    end
end
_merge(::Tuple{}, t2::Tuple, lt, by, rev) = t2
_merge(t1::Tuple, ::Tuple{}, lt, by, rev) = t1
_merge(::Tuple{}, ::Tuple{}, lt, by, rev) = ()