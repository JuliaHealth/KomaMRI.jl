function get_spin_coords(
    motion::ArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    Ns = size(motion.dx)[1]
    itpx, itpy, itpz = get_itp_functions(motion, Ns)
    if x isa CuArray
        itpx = adapt(CuArray{Float32}, itpx)
        itpy = adapt(CuArray{Float32}, itpy)
        itpz = adapt(CuArray{Float32}, itpz)  # Problem: too many CPU -> GPU transfers
    end
    return (x, y, z) .+ get_itp_results(itpx, itpy, itpz, Ns)
end