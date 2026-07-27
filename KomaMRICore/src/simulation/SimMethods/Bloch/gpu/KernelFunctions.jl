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

@inline function reduce_coils_serial!(
    sig_output, sig_r, sig_i, sens, sig_group_r, sig_group_i,
    i, i_l, i_g, s_idx, ADC_idx, N_spins, N_coils, N_adc, N, T, active,
    ::Val{MOTION}, ::Val{HAS_SENS}, warp_reduction,
) where {MOTION,HAS_SENS}
    coil = 1u32
    while coil <= N_coils
        coil_r, coil_i = sig_r, sig_i
        if active && HAS_SENS
            sens_idx = MOTION ? i + (s_idx - 1u32) * N_spins : i
            sens_r, sens_i = reim(sens[sens_idx, coil])
            coil_r, coil_i = (
                coil_r * sens_r - coil_i * sens_i,
                coil_r * sens_i + coil_i * sens_r,
            )
        end
        coil_r, coil_i = reduce_signal!(
            coil_r, coil_i, sig_group_r, sig_group_i, i_l, N, T, warp_reduction,
        )
        if i_l == 1u32
            sig_output[i_g, ADC_idx + (coil - 1u32) * N_adc] =
                complex(coil_r, coil_i)
        end
        coil += 1u32
    end
    return nothing
end

@inline function reduce_coils_tiled!(
    sig_output, sig_r, sig_i, sens, sig_group_r, sig_group_i,
    i, i_l, i_g, s_idx, ADC_idx, N_spins, N_coils, N_adc, N, T, active,
    ::Val{MOTION},
) where {MOTION}
    lane = (i_l - 1u32) % 32u32 + 1u32
    subgroup = (i_l - 1u32) ÷ 32u32 + 1u32
    nsubgroups = UInt32(N) ÷ 32u32
    first_coil = 1u32
    while first_coil <= N_coils
        tile_size = min(nsubgroups, N_coils - first_coil + 1u32)
        tile_coil = 1u32
        while tile_coil <= tile_size
            coil = first_coil + tile_coil - 1u32
            coil_r, coil_i = sig_r, sig_i
            if active
                sens_idx = MOTION ? i + (s_idx - 1u32) * N_spins : i
                sens_r, sens_i = reim(sens[sens_idx, coil])
                coil_r, coil_i = (
                    coil_r * sens_r - coil_i * sens_i,
                    coil_r * sens_i + coil_i * sens_r,
                )
            end
            coil_r, coil_i = reduce_warp(coil_r, coil_i)
            if lane == 1u32
                scratch_idx = (tile_coil - 1u32) * nsubgroups + subgroup
                sig_group_r[scratch_idx] = coil_r
                sig_group_i[scratch_idx] = coil_i
            end
            tile_coil += 1u32
        end
        @synchronize()

        if subgroup <= tile_size
            scratch_idx = (subgroup - 1u32) * nsubgroups + lane
            coil_r = lane <= nsubgroups ? sig_group_r[scratch_idx] : zero(T)
            coil_i = lane <= nsubgroups ? sig_group_i[scratch_idx] : zero(T)
            coil_r, coil_i = reduce_warp(coil_r, coil_i)
            if lane == 1u32
                coil = first_coil + subgroup - 1u32
                sig_output[i_g, ADC_idx + (coil - 1u32) * N_adc] =
                    complex(coil_r, coil_i)
            end
        end
        @synchronize()
        first_coil += nsubgroups
    end
    return nothing
end

## COV_EXCL_STOP
