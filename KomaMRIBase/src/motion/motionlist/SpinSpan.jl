abstract type AbstractSpinSpan end

# AllSpins
struct AllSpins <: AbstractSpinSpan end

Base.getindex(spins::AllSpins, p::AbstractVector) = Colon(), spins
Base.view(spins::AllSpins, p::AbstractVector) = Colon(), spins

Base.getindex(x, p::AllSpins) = x
Base.getindex(x, p::AllSpins, q) = x[:, q]
Base.view(x, p::AllSpins) = x
Base.view(x, p::AllSpins, q) = @view(x[:, q])

get_idx(spins::AllSpins) = Colon()
has_spins(spins::AllSpins) = true

# SpinRange
@with_kw struct SpinRange <: AbstractSpinSpan 
    range::AbstractVector
end

SpinRange(c::Colon) = AllSpins()
SpinRange(range::BitVector) = SpinRange(findall(x->x==true, range))

function Base.getindex(spins::SpinRange, p::AbstractVector)
    idx = get_idx(spins.range, p)
    spin_range = SpinRange(spins.range[idx])
    return idx, spin_range
end
function Base.view(spins::SpinRange, p::AbstractVector)
    idx = get_idx(spins.range, p)
    spin_range = SpinRange(@view(spins.range[idx]))
    return idx, spin_range
end

Base.getindex(spins::SpinRange, b::BitVector) = spins[findall(x->x==true, b)]
Base.view(spins::SpinRange, b::BitVector) = @view(spins[findall(x->x==true, b)])

Base.getindex(x, p::SpinRange) = x[p.range]
Base.getindex(x, p::SpinRange, q) = x[p.range, q]
Base.view(x, p::SpinRange) = @view(x[p.range])
Base.view(x, p::SpinRange, q) = @view(x[p.range, q])

Base.:(==)(sr1::SpinRange, sr2::SpinRange) = sr1.range == sr2.range

get_idx(spins::SpinRange) = spins.range
has_spins(spins::SpinRange) = length(spins.range) > 0

# Auxiliary functions
function get_idx(spin_range::AbstractVector, p::AbstractVector)
    idx = findall(x -> x in p, spin_range)
    return (length(idx) > 0 && idx == collect(first(idx):last(idx))) ? (first(idx):last(idx)) : idx
end

expand(sr::SpinRange, Ns::Int) = sr
expand(sr::AllSpins, Ns::Int) = SpinRange(1:Ns)