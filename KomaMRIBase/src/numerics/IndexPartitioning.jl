# ===========================================================================
# Numerical Index Partitioning
# ===========================================================================

"""
    array_of_ranges = kfoldperm(N, k; breaks=[])

Divides a list of indices from 1 to `N` into `k` groups.

# Arguments
- `N`: (`::Integer`) number of elements to be ordered
- `k`: (`::Integer`) number of groups to divide the `N` elements.

# Keywords
- `breaks`: (`::Vector{<:Integer}`, `=[]`) array of indices where predefined breakpoints are
    placed.

# Returns
- `array_of_ranges`: (`::Vector{UnitRange{<:Integer}}`) array containing ranges of different
    groups. The target is `k` groups, but this could increase by adding elements to the
    `breaks` input array
"""
function kfoldperm(N, k; breaks=[])
	k = min(N, k)
	n, r = divrem(N, k) #N >= k, N < k
	b = collect(1:n:N+1)
	Nb = length(b)
	for i in 1:Nb
		b[i] += i > r ? r : i-1
	end
	b = sort(unique(append!(b, breaks)))
	Nbreaks = length(b) - Nb
	p = 1:N
	return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k+Nbreaks]] #TODO: use RF starts and ends differently to remove PATCH in run_sim_time_iter
end
