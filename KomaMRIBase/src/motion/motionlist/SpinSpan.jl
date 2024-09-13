abstract type AbstractSpinSpan end

"""
    allspins = AllSpins()

AllSpins struct. It is a specialized type that inherits from AbstractSpinSpan
and is used to cover all the spins of a phantom.

# Returns
- `allspins`: (`::AllSpins`) AllSpins struct

# Examples
```julia-repl
julia> allspins = AllSpins()
```
"""
struct AllSpins <: AbstractSpinSpan end

# Functions
Base.getindex(spins::AllSpins, p) = p, spins
Base.view(spins::AllSpins, p) = p, spins
get_idx(spins::AllSpins) = Colon()
expand(sr::AllSpins, Ns::Int) = SpinRange(1:Ns)

"""
    spinrange = SpinRange(range)

SpinRange struct. It is a specialized type that inherits from AbstractSpinSpan
and is used to select a certain range and number of spins.

# Arguments
- `range`: (`::AbstractVector`) spin id's. This argument can be a Range, a Vector or a BitVector

# Returns
- `spinrange`: (`::SpinRange`) SpinRange struct

# Examples
```julia-repl
julia> spinrange = SpinRange(1:10)
julia> spinrange = SpinRange([1, 3, 5, 7])
julia> spinrange = SpinRange(obj.x .> 0)
```
"""
@with_kw struct SpinRange <: AbstractSpinSpan 
    range::AbstractRange
end

# Constructors
SpinRange(c::Colon) = AllSpins()
SpinRange(b::BitVector) = SpinRange(findall(x->x==true, b))
SpinRange(v::Vector) = begin
    step = v[2] - v[1]
    r = step == 1 ? (v[1]:v[end]) : (v[1]:step:v[end])
    @assert r == v "Cannot create a SpinRange with indices that are not evenly spaced (e.g., [1, 3, 4])."
    return SpinRange(r)
end

# Functions
function Base.getindex(spins::SpinRange, p)
    idx = intersect_idx(spins.range, p)
    return intersect_idx(p, spins.range), SpinRange(idx)
end
function Base.view(spins::SpinRange, p)
    idx = intersect_idx(spins.range, p)
    return intersect_idx(p, spins.range), SpinRange(idx)
end
Base.:(==)(sr1::SpinRange, sr2::SpinRange) = sr1.range == sr2.range
Base.length(sr::SpinRange) = length(sr.range)
get_idx(spins::SpinRange) = spins.range
expand(sr::SpinRange, Ns::Int) = sr
intersect_idx(a, b) = findall(x -> x in a, b)
