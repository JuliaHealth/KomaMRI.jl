## COV_EXCL_START

@kernel unsafe_indices=true inbounds=true function precession_kernel!(
    sig_output::AbstractMatrix{Complex{T}}, 
    M_xy, M_z, 
    @Const(p_x), @Const(p_y), @Const(p_z), @Const(p_ΔBz), @Const(p_T1), @Const(p_T2), @Const(p_ρ), N_spins,
    @Const(s_Gx), @Const(s_Gy), @Const(s_Gz), @Const(s_Δt), @Const(s_ADC), s_length,
    ::Val{MOTION}, ::Val{USE_WARP_REDUCTION},
) where {T, MOTION, USE_WARP_REDUCTION}

    @uniform N = @groupsize()[1]
    i_l = @index(Local, Linear)
    i_g = @index(Group, Linear)
    i = (i_g - 1u32) * UInt32(N) + i_l

    sig_group_r = @localmem T USE_WARP_REDUCTION ? 32 : N
    sig_group_i = @localmem T USE_WARP_REDUCTION ? 32 : N
    sig_r = zero(T)
    sig_i = zero(T)
    
    Mxy_r = zero(T)
    Mxy_i = zero(T)
    t = zero(T)
    ϕ = zero(T)
    ΔBz = zero(T)
    T2 = zero(T)
    x = zero(T)
    y = zero(T)
    z = zero(T)
    Bz_prev = zero(T)
    Bz_next = zero(T)

    if i <= N_spins
        Mxy_r, Mxy_i = reim(M_xy[i])
        sig_r = Mxy_r
        sig_i = Mxy_i
        ΔBz = p_ΔBz[i]
        T2 = p_T2[i]
        x, y, z = get_spin_coordinates(p_x, p_y, p_z, i, 1)
        Bz_prev = x * s_Gx[1] + y * s_Gy[1] + z * s_Gz[1] + ΔBz
    end

    ADC_idx = 1u32
    if s_ADC[1]
        sig_r, sig_i = reduce_signal!(sig_r, sig_i, sig_group_r, sig_group_i, i_l, N, T, Val(USE_WARP_REDUCTION))
        if i_l == 1u32
            sig_output[i_g, 1] = complex(sig_r, sig_i)
        end
        ADC_idx += 1u32
    end

    s_idx = 2u32
    while s_idx <= s_length
        if i <= N_spins
            if MOTION
                x, y, z = get_spin_coordinates(p_x, p_y, p_z, i, s_idx)
            end

            Δt = s_Δt[s_idx-1]
            t += Δt
            Bz_next = x * s_Gx[s_idx] + y * s_Gy[s_idx] + z * s_Gz[s_idx] + ΔBz
            ϕ += (Bz_prev + Bz_next) * T(-π * γ) * Δt
        end

        if s_idx < s_length && s_ADC[s_idx]
            if i <= N_spins
                ΔT2 = exp(-t / T2)
                cis_ϕ_i = sin(ϕ)
                cis_ϕ_r = cos(ϕ)
                sig_r = ΔT2 * (Mxy_r * cis_ϕ_r - Mxy_i * cis_ϕ_i)
                sig_i = ΔT2 * (Mxy_r * cis_ϕ_i + Mxy_i * cis_ϕ_r)
            end
            sig_r, sig_i = reduce_signal!(sig_r, sig_i, sig_group_r, sig_group_i, i_l, N, T, Val(USE_WARP_REDUCTION))
            if i_l == 1u32
                sig_output[i_g, ADC_idx] = complex(sig_r, sig_i)
            end
            ADC_idx += 1u32
        end
        
        Bz_prev = Bz_next
        s_idx += 1u32
    end

    if i <= N_spins
        ΔT1 = exp(-t / p_T1[i])
        ΔT2 = exp(-t / T2)
        cis_ϕ_i = sin(ϕ)
        cis_ϕ_r = cos(ϕ)
        M_xy[i] = complex(ΔT2 * (Mxy_r * cis_ϕ_r - Mxy_i * cis_ϕ_i), ΔT2 * (Mxy_r * cis_ϕ_i + Mxy_i * cis_ϕ_r))
        M_z[i] = M_z[i] * ΔT1 + p_ρ[i] * (T(1) - ΔT1)
    end
end

## COV_EXCL_STOP