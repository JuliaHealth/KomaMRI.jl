function index_samples(x, y)
    valid = Int32[]
    gaps = Int32[]
    for i in eachindex(x)
        if !ismissing(x[i]) && !ismissing(y[i]) && isfinite(x[i]) && isfinite(y[i])
            push!(valid, i)
        else
            push!(gaps, i)
        end
    end
    issorted(view(x, valid)) || sort!(valid; by=Base.Fix1(getindex, x))
    return (; valid, gaps)
end

function selected_samples(x, y, samples, range, bins)
    times = view(x, samples.valid)
    isempty(times) && return (Int[], 0)
    first_inside = searchsortedfirst(times, range[1])
    last_inside = searchsortedlast(times, range[2])
    available = max(0, last_inside - first_inside + 1)
    first_point = max(1, first_inside - 1)
    last_point = min(length(times), last_inside + 1)
    selected = Int[]
    if last_point - first_point + 1 <= 4bins
        append!(selected, view(samples.valid, first_point:last_point))
    else
        edges = LinRange(range[1], range[2], bins + 1)
        for bin in 1:bins
            a = bin == 1 ? first_point : searchsortedfirst(times, edges[bin])
            b = bin == bins ? last_point : searchsortedfirst(times, edges[bin + 1]) - 1
            a > b && continue
            imin = imax = Int(samples.valid[a])
            @inbounds for j in a:b
                i = Int(samples.valid[j])
                y[i] < y[imin] && (imin = i)
                y[i] > y[imax] && (imax = i)
            end
            for i in (Int(samples.valid[a]), imin, imax, Int(samples.valid[b]))
                # Keep both sides of event boundaries, including vertical ADC window edges.
                gap = searchsortedlast(samples.gaps, i)
                event_start = gap == 0 ? 1 : samples.gaps[gap] + 1
                event_end = gap == length(samples.gaps) ? length(x) : samples.gaps[gap + 1] - 1
                append!(selected, max(event_start, samples.valid[first_point]):
                    min(event_start + 1, event_end, samples.valid[last_point]))
                append!(selected, max(event_start, event_end - 1, samples.valid[first_point]):
                    min(event_end, samples.valid[last_point]))
                for j in max(1, i - 1):min(length(x), i + 1)
                    !ismissing(x[j]) && !ismissing(y[j]) && isfinite(x[j]) && isfinite(y[j]) &&
                        push!(selected, j)
                end
            end
        end
    end
    sort!(unique!(selected))
    indices = Int[]
    previous = 0
    for i in selected
        if previous != 0 && searchsortedlast(samples.gaps, i) != searchsortedlast(samples.gaps, previous)
            push!(indices, 0)
        end
        push!(indices, i)
        previous = i
    end
    return indices, available
end


function summarize_blocks(context, channels, data, layout, config, block_starts, full_range, read_block; break_blocks)
    nblocks = length(block_starts) - 1
    counts = zeros(Int32, nblocks, length(channels))
    offsets = ones(Int, nblocks + 1, length(channels))
    fields = map(eachindex(channels)) do j
        prototype = read_block(context, channels[j], 1)
        names = Tuple(name for name in keys(prototype) if haskey(data[j].fields, name))
        columns = NamedTuple{names}(map(names) do name
            Union{Missing,eltype(getproperty(prototype, name))}[]
        end)
        for block in 1:nblocks
            samples = read_block(context, channels[j], block)
            n = length(samples.x)
            counts[block, j] = n
            if n > 0
                extremal = (1, argmin(samples.y), argmax(samples.y), n)
                indices = sort!(unique([i for k in extremal for i in max(1, k - 1):min(n, k + 1)]))
                for name in names
                    column = getproperty(columns, name)
                    append!(column, getproperty(samples, name)[indices])
                    break_blocks[j] && push!(column, missing)
                end
            end
            offsets[block + 1, j] = length(columns.x) + 1
        end
        for name in names
            data[j][name] = getproperty(columns, name)
        end
        names
    end
    return (; context, channels, data, layout, config, counts, offsets, fields, block_starts, full_range, read_block, break_blocks)
end

function block_window_data(source, range, width, visibility)
    bins = clamp(floor(Int, width - 70), 24, 1000)
    starts = source.block_starts
    last_block = length(starts) - 1
    # Only blocks crossing pixel-bin edges need expanded samples to locate the local extrema.
    edges = LinRange(range[1], range[2], bins + 1)
    expanded = Set(clamp(searchsortedlast(starts, t), 1, last_block) for t in edges)
    lo = clamp(searchsortedlast(starts, range[1]), 1, last_block)
    hi = clamp(searchsortedlast(starts, range[2]), 1, last_block)
    available = 0
    shown = 0
    data = [copy(trace.fields) for trace in source.data]
    for j in eachindex(source.data)
        template = source.data[j]
        output = data[j]
        fields = source.fields[j]
        candidates = NamedTuple{fields}(empty(template[field]) for field in fields)
        if visibility[j] === true
            active = findall(>(0), @view(source.counts[:, j]))
            first_active = max(1, searchsortedfirst(active, lo) - 1)
            last_active = min(length(active), searchsortedlast(active, hi) + 1)
            full_detail = sum(@view(source.counts[lo:hi, j])) <= 4bins
            for block in @view(active[first_active:last_active])
                if full_detail || block in expanded || block < lo || block > hi
                    samples = source.read_block(source.context, source.channels[j], block)
                    available += count(t -> range[1] <= t <= range[2], samples.x)
                    for field in fields
                        append!(getproperty(candidates, field), getproperty(samples, field))
                        source.break_blocks[j] && push!(getproperty(candidates, field), missing)
                    end
                else
                    available += source.counts[block, j]
                    indices = source.offsets[block, j]:(source.offsets[block + 1, j] - 1)
                    for field in fields
                        append!(getproperty(candidates, field), @view(template[field][indices]))
                    end
                end
            end
        end
        indices, _ = selected_samples(candidates.x, candidates.y,
            index_samples(candidates.x, candidates.y), range, bins)
        visibility[j] == "legendonly" && (indices = [0])
        for field in fields
            values = getproperty(candidates, field)
            output[field] = [i == 0 ? nothing : values[i] for i in indices]
        end
        output[:mode] = get(output, :mode, "lines") == "line" ? "lines" : get(output, :mode, "lines")
        output[:uid] = "trace-$(j)"
        output[:visible] = visibility[j]
        shown += count(i -> i != 0 && range[1] <= candidates.x[i] <= range[2], indices)
    end
    return (; data, source.layout, source.config, source.full_range, shown, available)
end
