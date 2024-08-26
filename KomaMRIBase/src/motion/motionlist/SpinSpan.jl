abstract type AbstractSpinSpan end

"""
    as = AllSpin()

AllSpin struct. (...)  
"""
struct AllSpins <: AbstractSpinSpan end

Base.getindex(spins::AllSpins, p::AbstractVector) = p, spins
Base.view(spins::AllSpins, p::AbstractVector) = p, spins

get_idx(spins::AllSpins) = Colon()
has_spins(spins::AllSpins) = true

"""
    sr = SpinRange(range)

SpinRange struct. (...)  
"""
@with_kw struct SpinRange <: AbstractSpinSpan 
    range::AbstractVector
end

SpinRange(c::Colon) = AllSpins()
SpinRange(range::BitVector) = SpinRange(findall(x->x==true, range))

function Base.getindex(spins::SpinRange, p::AbstractVector)
    idx = get_idx(spins.range, p)
    return get_idx(p, spins.range), SpinRange(idx)
end
function Base.view(spins::SpinRange, p::AbstractVector)
    idx = get_idx(spins.range, p)
    return get_idx(p, spins.range), SpinRange(idx)
end

Base.getindex(spins::SpinRange, b::BitVector) = spins[findall(x->x==true, b)]
Base.view(spins::SpinRange, b::BitVector) = @view(spins[findall(x->x==true, b)])

Base.:(==)(sr1::SpinRange, sr2::SpinRange) = sr1.range == sr2.range

Base.length(sr::SpinRange) = length(sr.range)

get_idx(spins::SpinRange) = spins.range
has_spins(spins::SpinRange) = length(spins.range) > 0

# Auxiliary functions
function get_idx(spin_range::AbstractVector, p::AbstractVector)
    idx = findall(x -> x in spin_range, p)
    return (length(idx) > 0 && idx == collect(first(idx):last(idx))) ? (first(idx):last(idx)) : idx
end

expand(sr::SpinRange, Ns::Int) = sr
expand(sr::AllSpins, Ns::Int) = SpinRange(1:Ns)