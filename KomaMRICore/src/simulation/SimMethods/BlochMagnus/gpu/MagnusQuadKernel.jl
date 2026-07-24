## COV_EXCL_START

@kernel unsafe_indices=true inbounds=true function excitation_kernel!(
    sig_output::AbstractMatrix{Complex{T}},
    M_xy::AbstractVector{Complex{T}}, M_z,
    @Const(p_x), @Const(p_y), @Const(p_z), @Const(p_ΔBz), @Const(p_T1), @Const(p_T2), @Const(p_ρ), N_spins,
    @Const(s_Gx), @Const(s_Gy), @Const(s_Gz), @Const(s_Δt), @Const(s_Δf), @Const(s_B1), @Const(s_ψ), @Const(s_ADC), s_length,
    ::Val{MOTION}, ::Val{USE_WARP_REDUCTION}, ::Val{HAS_ADC}, ::Val{SKIP_RELAXATION},
    sim_method::SM
) where {T, MOTION, USE_WARP_REDUCTION, HAS_ADC, SKIP_RELAXATION, SM<:Union{BlochMagnusQuad2,BlochMagnusQuad4}}

    @uniform N = @groupsize()[1]
    i_l = @index(Local, Linear)
    i_g = @index(Group, Linear)
    i = (i_g - 1u32) * UInt32(N) + i_l

    sig_group_r = @localmem T HAS_ADC ? (USE_WARP_REDUCTION ? 32 : N) : 1
    sig_group_i = @localmem T HAS_ADC ? (USE_WARP_REDUCTION ? 32 : N) : 1

    active = i <= N_spins
    Mxy_r = zero(T)
    Mxy_i = zero(T)
    Mz = zero(T)
    ρ = zero(T)
    ΔBz = zero(T)
    T1 = T(1)
    T2 = T(1)
    x0 = zero(T)
    y0 = zero(T)
    z0 = zero(T)
    Bx_0 = zero(T)
    By_0 = zero(T)
    Bz_0 = zero(T)

    if active
        ΔBz = p_ΔBz[i]
        Mxy_r, Mxy_i = reim(M_xy[i])
        Mz = M_z[i]
        if !SKIP_RELAXATION
            ρ = p_ρ[i]
            T1 = p_T1[i]
            T2 = p_T2[i]
        end

        ψ_start = s_ψ[1]
        if !iszero(ψ_start)
            sin_ψ, cos_ψ = sincos(ψ_start)
            Mxy_r, Mxy_i = Mxy_r * cos_ψ + Mxy_i * sin_ψ, Mxy_i * cos_ψ - Mxy_r * sin_ψ
        end

        x0, y0, z0 = get_spin_coordinates(p_x, p_y, p_z, i, 1)
        Bx_0, By_0 = reim(s_B1[1])
        Bz_0 = x0 * s_Gx[1] + y0 * s_Gy[1] + z0 * s_Gz[1] + ΔBz - s_Δf[1] / T(γ)
    end

    ADC_idx = 1u32
    s_idx = 1u32
    while s_idx + 1u32 < s_length
        s_mid = s_idx + 1u32
        s_end = s_mid + 1u32

        if active
            xm, ym, zm = MOTION ? get_spin_coordinates(p_x, p_y, p_z, i, s_mid) : (x0, y0, z0)
            x1, y1, z1 = MOTION ? get_spin_coordinates(p_x, p_y, p_z, i, s_end) : (x0, y0, z0)

            Bx_m, By_m = reim(s_B1[s_mid])
            Bz_m = xm * s_Gx[s_mid] + ym * s_Gy[s_mid] + zm * s_Gz[s_mid] + ΔBz - s_Δf[s_mid] / T(γ)
            Bx_1, By_1 = reim(s_B1[s_end])
            Bz_1 = x1 * s_Gx[s_end] + y1 * s_Gy[s_end] + z1 * s_Gz[s_end] + ΔBz - s_Δf[s_end] / T(γ)

            Δt = s_Δt[s_idx] + s_Δt[s_mid]
            θx, θy, θz = rotation_vector(
                Bx_0, By_0, Bz_0,
                Bx_m, By_m, Bz_m,
                Bx_1, By_1, Bz_1,
                Δt,
                sim_method,
            )
            M_norm = mag_norm(T, Mxy_r, Mxy_i, Mz)
            Mxy_r, Mxy_i, Mz = rotate_magnetization(θx, θy, θz, Mxy_r, Mxy_i, Mz, T)
            Mxy_r, Mxy_i, Mz = restore_mag_norm(M_norm, Mxy_r, Mxy_i, Mz) # For reduced float precision only.

            if !SKIP_RELAXATION
                E1 = exp(-Δt / T1)
                E2 = exp(-Δt / T2)
                Mxy_r *= E2
                Mxy_i *= E2
                Mz = Mz * E1 + ρ * (T(1) - E1)
            end

            x0, y0, z0 = x1, y1, z1
            Bx_0, By_0, Bz_0 = Bx_1, By_1, Bz_1
        end

        if HAS_ADC && s_ADC[s_end]
            sig_r, sig_i = reduce_signal!(Mxy_r, Mxy_i, sig_group_r, sig_group_i, i_l, N, T, Val(USE_WARP_REDUCTION))
            if i_l == 1u32
                sig_output[i_g, ADC_idx] = complex(sig_r, sig_i)
            end
            ADC_idx += 1u32
        end

        s_idx = s_end
    end

    if active
        ψ_end = s_ψ[s_length]
        if !iszero(ψ_end)
            sin_ψ, cos_ψ = sincos(ψ_end)
            Mxy_r, Mxy_i = Mxy_r * cos_ψ - Mxy_i * sin_ψ, Mxy_r * sin_ψ + Mxy_i * cos_ψ
        end
        M_xy[i] = complex(Mxy_r, Mxy_i)
        M_z[i] = Mz
    end
end

## COV_EXCL_STOP
