using KernelAbstractions: @kernel, @Const, @index, @uniform, @localmem, @synchronize, @groupsize, @subgroupreduce
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

@inline reduce_subgroup(val_r, val_i) =
    reim(@subgroupreduce(+, complex(val_r, val_i), zero(complex(val_r, val_i))))

# Reduce every sub-group with shuffles, then the per-sub-group sums in every sub-group.
@inline function reduce_signal!(sig_r, sig_i, sig_group_r, sig_group_i, i_l, N, T, ::Val{true})
    sig_r, sig_i = reduce_subgroup(sig_r, sig_i)

    subgroup = KI.get_sub_group_id(UInt32)
    lane = KI.get_sub_group_local_id(UInt32)
    if lane == 1u32
        @inbounds sig_group_r[subgroup] = sig_r
        @inbounds sig_group_i[subgroup] = sig_i
    end

    @synchronize()

    sig_r = zero(T)
    sig_i = zero(T)
    i = lane
    while i <= KI.get_num_sub_groups(UInt32)
        @inbounds sig_r += sig_group_r[i]
        @inbounds sig_i += sig_group_i[i]
        i += KI.get_sub_group_size(UInt32)
    end
    sig_r, sig_i = reduce_subgroup(sig_r, sig_i)

    # All sub-groups must finish reading before the next ADC overwrites scratch.
    @synchronize()
    return sig_r, sig_i
end

# GPU-kernel sensitivity lookup: matrices are precomputed maps, while receiver
# models evaluate one spin and coil directly at the current coordinates.
@inline function KomaMRIBase.get_sens(
    sens::AbstractMatrix,
    _positions::Tuple{AbstractVector,AbstractVector,AbstractVector},
    spin, _s_idx, _N_spins, coil,
)
    return @inbounds sens[spin, coil]
end

@inline function KomaMRIBase.get_sens(
    sens::AbstractMatrix,
    _positions::Tuple{AbstractMatrix,AbstractMatrix,AbstractMatrix},
    spin, s_idx, N_spins, coil,
)
    return @inbounds sens[spin + (s_idx - 1u32) * N_spins, coil]
end

@inline function KomaMRIBase.get_sens(
    receiver, positions, spin, s_idx, _N_spins, coil,
)
    x, y, z = positions
    position = get_spin_coordinates(x, y, z, spin, s_idx)
    return get_sens(receiver, position, coil)
end

# Fallback for backends without sub-group shuffles.
@inline function reduce_signal_per_coil!(
    sig_output, sig_r, sig_i, receiver, sig_group_r, sig_group_i,
    positions, s_idx,
    i_l, i_g, ADC_idx, N_spins, N_coils, N_adc, N, T,
    ::Val{false},
)
    spin = (i_g - 1u32) * UInt32(N) + i_l
    active = spin <= N_spins
    coil = 1u32
    while coil <= N_coils
        coil_r, coil_i = sig_r, sig_i
        if active
            sens_r, sens_i = reim(get_sens(
                receiver, positions, spin, s_idx, N_spins, coil,
            ))
            coil_r, coil_i = (
                coil_r * sens_r - coil_i * sens_i,
                coil_r * sens_i + coil_i * sens_r,
            )
        end
        coil_r, coil_i = reduce_signal!(
            coil_r, coil_i, sig_group_r, sig_group_i, i_l, N, T, Val(false),
        )
        if i_l == 1u32
            @inbounds sig_output[i_g, ADC_idx + (coil - 1u32) * N_adc] =
                complex(coil_r, coil_i)
        end
        coil += 1u32
    end
    return nothing
end

# Assign one coil to each sub-group when sub-group shuffles are available.
@inline function reduce_signal_per_coil!(
    sig_output, sig_r, sig_i, receiver, sig_group_r, sig_group_i,
    positions, s_idx,
    i_l, i_g, ADC_idx, N_spins, N_coils, N_adc, N, T,
    ::Val{true},
)
    @inbounds sig_group_r[i_l] = sig_r
    @inbounds sig_group_i[i_l] = sig_i
    @synchronize()

    lane = KI.get_sub_group_local_id(UInt32)
    subgroup = KI.get_sub_group_id(UInt32)
    nsubgroups = KI.get_num_sub_groups(UInt32)
    width = KI.get_sub_group_size(UInt32)

    # Every sub-group runs the same number of rounds: some backends (PoCL) need sub-group
    # operations in work-group-uniform control flow.
    first_coil = 1u32
    while first_coil <= N_coils
        coil = first_coil + subgroup - 1u32
        coil_r = zero(T)
        coil_i = zero(T)
        local_spin = lane
        while coil <= N_coils && local_spin <= N
            spin = (i_g - 1u32) * UInt32(N) + local_spin
            if spin <= N_spins
                sens_r, sens_i = reim(get_sens(
                    receiver, positions, spin, s_idx, N_spins, coil,
                ))
                signal_r = @inbounds sig_group_r[local_spin]
                signal_i = @inbounds sig_group_i[local_spin]
                coil_r += signal_r * sens_r - signal_i * sens_i
                coil_i += signal_r * sens_i + signal_i * sens_r
            end
            local_spin += width
        end
        coil_r, coil_i = reduce_subgroup(coil_r, coil_i)
        if lane == 1u32 && coil <= N_coils
            @inbounds sig_output[i_g, ADC_idx + (coil - 1u32) * N_adc] =
                complex(coil_r, coil_i)
        end
        first_coil += nsubgroups
    end

    # All subgroups must finish reading before the next ADC overwrites scratch.
    @synchronize()
    return nothing
end

## COV_EXCL_STOP
