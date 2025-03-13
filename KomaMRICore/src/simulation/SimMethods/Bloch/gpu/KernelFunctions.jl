using KernelAbstractions: @kernel, @Const, @index, @uniform, @localmem, @synchronize, @groupsize
using KernelAbstractions.Extras: @unroll

## COV_EXCL_START

#Used for getting spin coordinates inside precession and excitation kernels
@inline function get_spin_coordinates(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, i::Integer, t::Integer) where {T<:Real} 
    @inbounds (x[i], y[i], z[i]) 
end
@inline function get_spin_coordinates(x::AbstractMatrix{T}, y::AbstractMatrix{T}, z::AbstractMatrix{T}, i::Integer, t::Integer) where {T<:Real} 
    @inbounds (x[i, t], y[i, t], z[i, t]) 
end

# Returns the next least power of two starting from n, used to calculate remaining indexes in the first step of a threadgroup-level reduction.
@inline function next_least_power_of_two(n)
    return  n < 2 ? 1 :
            n < 4 ? 2 :
            n < 8 ? 4 :
            n < 16 ? 8 :
            n < 32 ? 16 :
            n < 64 ? 32 :
            n < 128 ? 64 :
            n < 256 ? 128 :
            n < 512 ? 256 :
            n < 1024 ? 512 :
            1024
end

@inline function num_reduction_iterations(n)
    return  n == 2 ? 1 :
            n == 4 ? 2 :
            n == 8 ? 3 :
            n == 16 ? 4 :
            n == 32 ? 5 :
            n == 64 ? 6 :
            n == 128 ? 7 :
            n == 256 ? 8 :
            n == 512 ? 9 :
            10
end

@inline function reduce_signal!(sig_r, sig_i, i_l, N)
    N_closest = next_least_power_of_two(N)
    if N != N_closest
        R = UInt32(N - N_closest)
        if i_l <= R
            @inbounds sig_r[i_l] = sig_r[i_l] + sig_r[i_l + N_closest]
            @inbounds sig_i[i_l] = sig_i[i_l] + sig_i[i_l + N_closest]
        end
        @synchronize()
    end

    @unroll for k=1:num_reduction_iterations(N_closest)
        offset = N_closest >> k
        if i_l <= offset
            @inbounds sig_r[i_l] = sig_r[i_l] + sig_r[i_l + offset]
            @inbounds sig_i[i_l] = sig_i[i_l] + sig_i[i_l + offset]
        end
        @synchronize()
    end
end

## COV_EXCL_STOP