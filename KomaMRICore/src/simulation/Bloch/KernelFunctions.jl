using KernelAbstractions: @index, @kernel

"""
    cumsum2_kernel

Simple kernel function, computes the cumulative sum of each row of a matrix. Operates
in-place on the input matrix without allocating additional memory.

# Arguments
- 'A': matrix to compute cumsum on
"""
@kernel function cumsum_matrix_rows_kernel!(A)
    i = @index(Global)

    for k ∈ 2:size(A, 2)
        @inbounds A[i, k] += A[i, k-1]
    end
end

"""
    cumtrapz

A more efficient GPU implementation of the cumtrapz method defined in TrapezoidalIntegration.jl.
Uses a kernel to compute cumsum along the second dimension.

# Arguments
- `Δt`: (`1 x NΔt ::Matrix{Float64}`, `[s]`) delta time 1-row array
- `x`: (`Ns x (NΔt+1) ::Matrix{Float64}`, `[T]`) magnitude of the field Gx * x + Gy * y +
    Gz * z

# Returns
- `y`: (`Ns x NΔt ::Matrix{Float64}`, `[T*s]`) matrix where every column is the
    cumulative integration over time of (Gx * x + Gy * y + Gz * z) * Δt for every spin of a
    phantom
"""
function KomaMRIBase.cumtrapz(Δt::AbstractArray{T}, x::AbstractArray{T}, backend::KA.GPU) where {T<:Real}
    y = (x[:, 2:end] .+ x[:, 1:end-1]) .* (Δt / 2)
    cumsum_matrix_rows_kernel!(backend)(y, ndrange=size(y,1))
    KA.synchronize(backend)
    return y
end

#If on CPU, forward call to cumtrapz in KomaMRIBase
KomaMRIBase.cumtrapz(Δt::AbstractArray{T}, x::AbstractArray{T}, backend::KA.CPU) where {T<:Real} = KomaMRIBase.cumtrapz(Δt, x)