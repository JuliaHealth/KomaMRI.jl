function get_spin_coords(
    motion::ArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    Ns = size(motion.dx)[1]
    itpx, itpy, itpz = KomaMRIBase.get_itp_functions(motion, Ns)
    # To GPU
    itpx = adapt(CuArray{T}, itpx)
    itpy = adapt(CuArray{T}, itpy)
    itpz = adapt(CuArray{T}, itpz)  # Problem: too many CPU -> GPU transfers
    ux, uy, uz = KomaMRIBase.get_itp_results(itpx, itpy, itpz, t, Ns)
    return x .+ ux, y .+ uy, z .+ uz
end

