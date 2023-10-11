"""
For interpolating and extrapolating values
td: time vector of the the event defined by the user, it must be a non-decreasing time vector
ad: amplitude vector of the event defined by the user
Can be more than one amplitude at the same time
ts: sampler time vector, it is meant to be an non-decreasing time vector
"""
function interpolate(td::Vector{Float64}, ad::Vector{<:Number}, ts::Vector{Float64})

    # Copy the elementes of the time-event-values
    te = copy(td)
    ae = copy(ad)

    # Get the limits of the time-event-values region
    te0, ae0 = te[1], ae[1]
    tef, aef = te[end], ae[end]

    # Empty Vectors to be filled
    t, a = Float64[], Number[]  # Be careful with type stability here
    Nt = 0

    # Iterate over the sampling vector
    Nsi, Nei = 1, 0
    ion0, ionf = 0, 0
    isonfirst, isonlast = true, true
    for tsi in ts
        # Fill for time-samples previous to the time-event-values
        if tsi < te0
            append!(t, tsi); append!(a, ae0)
            Nt += 1
        # Fill for After the time-event-values region is covered
        elseif tef < tsi
            # Assign the index of the last index when is the event is on
            if isonlast
                ionf, isonlast = Nt, false
            end
            append!(t, tsi); append!(a, aef)
            Nt += 1
        # Fill for time-samples in the time-event-values region
        else
            # Initial condition for starting to consider the event values
            if isonfirst
                # Assign the index of the first index when is the event is on
                ion0, isonfirst = Nt + 1, false
                while te[1] < tsi
                    te0, ae0 = te[1], ae[1]
                    popfirst!(te); popfirst!(ae)
                end
            end
            # Fill with the exact values until all the samples at the same time are considered
            if !isempty(te) && tsi == te[1]
                Nsi, Nei = 1, 0
                while !isempty(te) && tsi == te[1]
                    te0, ae0 = te[1], ae[1]
                    #println("$tsi, $te0, $te")
                    append!(t, te0); append!(a, ae0)
                    popfirst!(te); popfirst!(ae)
                    Nei += 1
                    Nt += 1
                end
            # Fill with additional samples when more samples are needed at the same time
            elseif tsi == te0
                Nsi += 1
                #println("$Nsi, $Nei")
                if Nei < Nsi
                    append!(t, tsi); append!(a, ae0)
                    Nt += 1
                end
            # Fill with the linear interpolated value
            else
                ai = ae0 + (tsi-te0)*(ae[1]-ae0)/(te[1]-te0)
                #println("$tsi, $te0, $te")
                append!(t, tsi); append!(a, ai)
                Nt += 1
            end
        end
    end

    # Assign the index of the last index when is the event is on if it is not covered
    if !isonfirst && ionf == 0
        ionf = Nt
    end

    # Return the interpolation
    return (t = t, a = a, ion = (ion0, ionf))
end

"""
Returns a vector with the same time-samples from two non-decreasing time-vectors,
but it keeps the maximum number of time-samples found in both vectors when they are repeated.
The merged time-samples are in increasing order
"""
function mergetimes(t1::Vector{Float64}, t2::Vector{Float64})
    # Empty vector to be filled
    t = Float64[]
    # Define lengths and initial counters
    Nt1, Nt2 = length(t1), length(t2)
    i, j = 1, 1
    # Core iteration
    while i <= Nt1 && j <= Nt2
        if t1[i] < t2[j]
            push!(t, t1[i]); i += 1
        elseif t1[i] > t2[j]
            push!(t, t2[j]); j += 1
        else
            push!(t, t1[i]); i += 1; j += 1
        end
    end
    # Append any remaining elements (just one "while-loop" can occur)
    while i <= Nt1
        push!(t, t1[i]); i += 1
    end
    while j <= Nt2
        push!(t, t2[j]); j += 1
    end
    # Return the merged times
    return t
end

"""
Same as merge times, but it considers a "vector" of non-decreasing time-vectors
"""
function mergetimes(vtimes::Vector{Vector{Float64}})
    t = vtimes[1]
    for τ in vtimes[2:end]
        t = mergetimes(t, τ)
    end
    return t
end

function is_intersection(a::Tuple{Int64, Int64}, b::Tuple{Int64, Int64})
    return !(a[2] < b[1] || b[2] < a[1])
end

function interval_minmax(a::Tuple{Int64, Int64}, b::Tuple{Int64, Int64})
    return (min(a[1], b[1]), max(a[2], b[2]))
end

"""
It returns a vector of tuples
The returned vectors are in increasing order
"""
function interval_union(a::Tuple{Int64, Int64}, b::Tuple{Int64, Int64})
    if is_intersection(a, b)
        return [interval_minmax(a, b)]
    end
    return (a[2] < b[1]) ? [a; b] : [b; a]
end
function interval_union(a::Tuple{Int64, Int64}, b::Tuple{Int64, Int64}, c::Tuple{Int64, Int64})
    if is_intersection(a, b)
        #ab
        ab = interval_minmax(a, b)
        return interval_union(ab, c)
    else
        #a,b
        if is_intersection(a, c)
            #ac
            ac = interval_minmax(a, c)
            return interval_union(ac, b)
        else
            #a,c
            if is_intersection(b, c)
                #bc
                bc = interval_minmax(b, c)
                return interval_union(bc, a)
            else
                #b,c
                v = sort([a[1], b[1], c[1]])
                if v[1] == a[1]
                    return (v[2] == b[1]) ? [a; b; c] : [a; c; b]
                elseif v[1] == b[1]
                    return (v[2] == a[1]) ? [b; a; c] : [b; c; a]
                else
                    return (v[2] == a[1]) ? [c; a; b] : [c; b; a]
                end
            end
        end
    end
end

"""
It returns a vector of tuples
The returned vectors are in increasing order
"""
function interval_intersection(a::Tuple{Int64, Int64}, b::Tuple{Int64, Int64})
    if a[2] <= b[1]
        return [(a[1], a[2], 0); (b[1], b[2], 1)]
    elseif b[2] <= a[1]
        return [(b[1], b[2], 1); (a[1], a[2], 0)]
    elseif (a[1] < b[1]) && (a[2] < b[2])
        return [(a[1], b[1], 0); (b[1], a[2], 2); (a[2], b[2], 1)]
    elseif (b[1] < a[1]) && (b[2] < a[2])
        return [(b[1], a[1], 1); (a[1], b[2], 2); (b[2], a[2], 0)]
    elseif a[1] == b[1]
        if a[2] == b[2]
            return [(a[1], a[2], 2)]
        elseif a[2] < b[2]
            return [(a[1], a[2], 2); (a[2], b[2], 1)]
        elseif b[2] < a[2]
            return [(b[1], b[2], 2); (b[2], a[2], 0)]
        end
    elseif a[2] == b[2]
        if a[1] == b[1]
            return [(b[1], b[2], 2)]
        elseif b[1] < a[1]
            return [(b[1], a[1], 1); (a[1], a[2], 2)]
        elseif a[1] < b[1]
            return [(a[1], b[1], 0); (b[1], b[2], 2)]
        end
    elseif (b[1] < a[1]) && (a[2] < b[2])
        return [(b[1], a[1], 1); (a[1], a[2], 2); (a[2], b[2], 1)]
    elseif (a[1] < b[1]) && (b[2] < a[2])
        return [(a[1], b[1], 0); (b[1], b[2], 2); (b[2], a[2], 0)]
    end
end

"""
This computation asummes that Δtrf < Δtgr
"""
function rfgr_intersection(rfion::Tuple{Int64, Int64}, grion::Vector{Tuple{Int64, Int64}})
    Ngrion = length(grion)
    if Ngrion == 0
        return [(rfion[1], rfion[2], 1)]
    elseif Ngrion == 1
        return interval_intersection(grion[1], rfion)
    elseif Ngrion == 2
        if rfion[2] < grion[2][1]
            return [interval_intersection(grion[1], rfion); (grion[2][1], grion[2][2], 0)]
        else
            return interval_intersection(interval_union(grion[1], grion[2]), rfion)
        end
    else
        if rfion[2] < grion[3][1]
            return [interval_intersection(interval_union(grion[1], grion[2]), rfion), (grion[3][1], grion[3][2], 0)]
        else
            return interval_intersection(interval_union(grion[1], grion[2], grion[3]), rfion)
        end
    end
end

"""
For refining times, since it uses bion, it depends on rfgr_intersection(), so it assumes that Δtrf < Δtgr
"""
function refinetimes(tsc::Vector{Float64}, bion::Vector{Tuple{Int64, Int64, Int64}}, Δtgr::Float64, Δtrf::Float64)
    # If the on-intervals is empty, then no refinement is necessary
    if isempty(bion)
       return tsc
    end
    # Perform the refinement
    ts = Float64[]
    i, Ni = 1, length(tsc)
    j, Nj = 1, length(bion)
    Δtr = (bion[1][3] == 0) ? Δtgr : Δtrf
    while j <= Nj
        # Fill with critical times
        if i <= bion[j][1]
            push!(ts, tsc[i])
        # Fill with refined times if necesary
        elseif i < bion[j][2]
            (tsc[i] - tsc[i-1] > Δtr) ? append!(ts, Vector(reverse(tsc[i]:-Δtr:tsc[i-1]))) : push!(ts, tsc[i])
        # Last interval with refined times if necesary
        else
            (tsc[i] - tsc[i-1] > Δtr) ? append!(ts, Vector(reverse(tsc[i]:-Δtr:tsc[i-1]))) : push!(ts, tsc[i])
            j += 1  # Update the on-interval counter
            (j <= Nj) && (Δtr = (bion[j][3] == 0) ? Δtgr : Δtrf)
        end
        i += 1  # Update the time counter
    end
    # Fill the remaining samples
    while i <= Ni
        push!(ts, tsc[i])
        i += 1
    end
    return ts
end

"""
all the elements of tadc must be contained in ts
"""
function indices_adcon(ts::Vector{Float64}, tadc::Vector{Float64})
    iadcon = Int64[]
    i, j, Ni, Nj = 1, 1, length(ts), length(tadc)
    while j <= Nj
        if ts[i] == tadc[j]
            while i <= Ni && ts[i] == tadc[j]
                push!(iadcon, i)
                i += 1
            end
            j += 1
        else
            i += 1
        end
    end
    return iadcon
end
"""
all the elements of tadc must be contained in ts
"""
function mask_adcon(ts::Vector{Float64}, tadc::Vector{Float64})
    i, j, Ni, Nj = 1, 1, length(ts), length(tadc)
    mask = fill(false, Ni)
    while j <= Nj
        if ts[i] == tadc[j]
            while i <= Ni && ts[i] == tadc[j]
                mask[i] = true
                i += 1
            end
            j += 1
        else
            i += 1
        end
    end
    return mask
end

"""
"""
function block_limits(addfirst::Bool, addlast::Bool, ΔT::Float64)
    (!addfirst && !addlast) && return Float64[]
    (!addfirst)             && return [ΔT]
    (!addlast)              && return [0.]
    return [0.; ΔT]
end


"""
Returns the first index value when tx appears in ts or is the next closest index
tx must be contained in the ts samples
"""
function index_rfx(ts::Vector{Float64}, tx::Float64)
    ix = 0
    if isnan(tx) || isempty(ts)
        return ix
    end
    ix = 1
    while true
        if ts[ix] >= tx
            return ix
        end
        ix += 1
    end
end


"""
all the elements of tc must be contained in ts
"""
function indices_rfon(ts::Vector{Float64}, tc::Vector{Float64})
    irfon = Int64[]
    if isempty(tc)
        return irfon
    end
    tmin, tmax = tc[1], tc[end]
    # Find min index
    i = 1
    while true
        if ts[i] == tmin
            push!(irfon, i)
            break
        end
        i += 1
    end
    # Find max index
    i = length(ts)
    while true
        if ts[i] == tmax
            push!(irfon, i)
            break
        end
        i -= 1
    end
    return irfon
end

function mask_rfon(ts::Vector{Float64}, tc::Vector{Float64})
    mask = fill(false, length(ts))
    if isempty(tc)
        return mask
    end
    imin, imax = indices_rfon(ts, tc)
    mask[imin:imax] .= true
    return mask
end

"""
It does same as get_sim_ranges() but ignores the Minimal quantity of number of Blocks
"""
function simrangesold(sq_irfon::Vector{Vector{Int64}}, Nt::Int64)
    ########################################################################################
    # Remove the last sample, this is due to DiscretizedSequence get by ranges, this shoudn't be done
    Nt -= 1
    ########################################################################################
    parts, excitation_bool = UnitRange{Int}[], Bool[]
    i0, Ni = 1, length(sq_irfon)
    if sq_irfon[1][1] != 1
        push!(parts, (1:(sq_irfon[1][1]-1))); push!(excitation_bool, false)
        i0 = sq_irfon[1][1]
    end
    for i in eachindex(sq_irfon)
        push!(parts, (sq_irfon[i][1]:sq_irfon[i][2])); push!(excitation_bool, true)
        if i < Ni
            push!(parts, ((sq_irfon[i][2]+1):(sq_irfon[i+1][1]-1))); push!(excitation_bool, false)
        end
    end
    if sq_irfon[end][2] != Nt
        push!(parts, ((sq_irfon[end][2]+1):Nt)); push!(excitation_bool, false)
    end
    return parts, excitation_bool
end

"""
Get the simulation ranges when RF is on or not
"""
function simranges(sq_irfon::Vector{Vector{Int64}}, Nt::Int64)
    parts, excitation_bool = UnitRange{Int}[], Bool[]
    Ni = length(sq_irfon)
    if sq_irfon[1][1] != 1
        push!(parts, (1:sq_irfon[1][1])); push!(excitation_bool, false)
    end
    for i in eachindex(sq_irfon)
        push!(parts, (sq_irfon[i][1]:sq_irfon[i][2])); push!(excitation_bool, true)
        if i < Ni
            push!(parts, (sq_irfon[i][2]:sq_irfon[i+1][1])); push!(excitation_bool, false)
        end
    end
    if sq_irfon[end][2] != Nt
        push!(parts, (sq_irfon[end][2]:Nt)); push!(excitation_bool, false)
    end
    return parts, excitation_bool
end
