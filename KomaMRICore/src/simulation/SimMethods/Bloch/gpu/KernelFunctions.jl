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

@inline function reduce_signal!(sig_r, sig_i, sig_group_r, sig_group_i, i_l, N, T, ::Val{false})
    @inbounds sig_group_r[i_l] = sig_r
    @inbounds sig_group_i[i_l] = sig_i
    @synchronize()

    N_closest = next_least_power_of_two(N)
    if N != N_closest
        R = UInt32(N - N_closest)
        if i_l <= R
            @inbounds sig_group_r[i_l] = sig_group_r[i_l] + sig_group_r[i_l + N_closest]
            @inbounds sig_group_i[i_l] = sig_group_i[i_l] + sig_group_i[i_l + N_closest]
        end
        @synchronize()
    end

    @unroll for k=1:num_reduction_iterations(N_closest)
        offset = N_closest >> k
        if i_l <= offset
            @inbounds sig_group_r[i_l] = sig_group_r[i_l] + sig_group_r[i_l + offset]
            @inbounds sig_group_i[i_l] = sig_group_i[i_l] + sig_group_i[i_l + offset]
        end
        @synchronize()
    end

    return sig_group_r[i_l], sig_group_i[i_l]
end

@inline function reduce_warp(val1, val2)
    @unroll for k=0:4
        val1 = val1 + shfl_down(val1, 1u32 << k)
        val2 = val2 + shfl_down(val2, 1u32 << k)
    end
    return val1, val2
end

@inline function reduce_signal!(sig_r, sig_i, sig_group_r, sig_group_i, i_l, N, T, ::Val{true})
    sig_r, sig_i = reduce_warp(sig_r, sig_i)

    if i_l % 32u32 == 1u32
        @inbounds sig_group_r[i_l ÷ 32u32 + 1u32] = sig_r
        @inbounds sig_group_i[i_l ÷ 32u32 + 1u32] = sig_i
    end

    @synchronize()

    @inbounds sig_r = (i_l <= UInt32(N) ÷ 32u32) ? sig_group_r[i_l] : zero(T)
    @inbounds sig_i = (i_l <= UInt32(N) ÷ 32u32) ? sig_group_i[i_l] : zero(T)
    
    return reduce_warp(sig_r, sig_i)
end

@inline function get_Bz!(seq::DiscreteHigherOrderSequence{T, Val(-1)}, x, y, z, s_idx::UInt32) where {T<:Real}
    return x * seq.G[1, s_idx] + y * seq.G[2, s_idx] + z * seq.G[3, s_idx]
end

@inline function get_Bz!(seq::DiscreteHigherOrderSequence{T, Val(0)}, x, y, z, s_idx::UInt32) where {T<:Real}
    return seq.G[1, s_idx]
end

@inline function get_Bz!(seq::DiscreteHigherOrderSequence{T, Val(1)}, x, y, z, s_idx::UInt32) where {T<:Real}
    return seq.G[1, s_idx] + x * seq.G[2, s_idx] + y * seq.G[3, s_idx] + z * seq.G[4, s_idx]
end

@inline function get_Bz!(seq::DiscreteHigherOrderSequence{T, Val(2)}, x, y, z, s_idx::UInt32) where {T<:Real}    
    return seq.G[1, s_idx] + 
           seq.G[2, s_idx] * x + 
           seq.G[3, s_idx] * y + 
           seq.G[4, s_idx] * z +
           seq.G[5, s_idx] * x * y + 
           seq.G[6, s_idx] * y * z + 
           seq.G[7, s_idx] * (2.0*z^2 - x^2 - y^2) +
           seq.G[8, s_idx] * x * z + 
           seq.G[9, s_idx] * (x^2 - y^2)
end

## COV_EXCL_STOP