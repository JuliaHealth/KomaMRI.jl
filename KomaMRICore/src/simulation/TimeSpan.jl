function unit_time(t::AbstractArray{T}, ts::TimeRange{T}) where {T<:Real}
    if ts.t_start == ts.t_end
        return (t .>= ts.t_start) .* oneunit(T)
    else
        tmp = max.((t .- ts.t_start) ./ (ts.t_end - ts.t_start), zero(T))
        t = min.(tmp, oneunit(T))
        KA.synchronize(KA.get_backend(t))
        return t
    end
end