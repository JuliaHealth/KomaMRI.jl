## COV_EXCL_START

@kernel unsafe_indices=true inbounds=true function excitation_kernel!(
    sig_output::AbstractMatrix{Complex{T}},
    M_xy::AbstractVector{Complex{T}}, M_z,
    @Const(p_x), @Const(p_y), @Const(p_z), @Const(p_ΔBz), @Const(p_T1), @Const(p_T2), @Const(p_ρ), N_spins,
    @Const(s_Gx), @Const(s_Gy), @Const(s_Gz), @Const(s_Δt), @Const(s_Δf), @Const(s_B1), @Const(s_ψ), @Const(s_ADC), s_length,
    ::Val{MOTION}, ::Val{USE_WARP_REDUCTION}, ::Val{HAS_ADC},
    sim_method::SM
) where {T, MOTION, USE_WARP_REDUCTION, HAS_ADC, SM<:Union{BlochMagnusGL2,BlochMagnusGL4}}

    @uniform N = @groupsize()[1]
    i_l = @index(Local, Linear)
    i_g = @index(Group, Linear)
    i = (i_g - 1u32) * UInt32(N) + i_l

    inv_γ = inv(T(γ))
    sig_group_r = @localmem T HAS_ADC ? (USE_WARP_REDUCTION === false ? N : 32) : 1
    sig_group_i = @localmem T HAS_ADC ? (USE_WARP_REDUCTION === false ? N : 32) : 1

    active = i <= N_spins
    Mxy_r = zero(T)
    Mxy_i = zero(T)
    Mz = zero(T)
    ρ = zero(T)
    ΔBz = zero(T)
    T1 = T(1)
    T2 = T(1)
    neg_inv_T1 = T(-1)
    neg_inv_T2 = T(-1)
    x = zero(T)
    y = zero(T)
    z = zero(T)

    if active
        ΔBz = p_ΔBz[i]
        Mxy_r, Mxy_i = reim(M_xy[i])
        Mz = M_z[i]
        ρ = p_ρ[i]
        T1 = p_T1[i]
        T2 = p_T2[i]
        neg_inv_T1 = -inv(T1)
        neg_inv_T2 = -inv(T2)

        ψ_start = s_ψ[1]
        if !iszero(ψ_start)
            sin_ψ, cos_ψ = sincos(ψ_start)
            Mxy_r, Mxy_i = Mxy_r * cos_ψ + Mxy_i * sin_ψ, Mxy_i * cos_ψ - Mxy_r * sin_ψ
        end

        x, y, z = get_spin_coordinates(p_x, p_y, p_z, i, 1)
    end

    ADC_idx = 1u32
    s_idx = 1u32
    while s_idx + 2u32 < s_length
        s_minus = s_idx + 1u32
        s_plus = s_minus + 1u32
        s_end = s_plus + 1u32

        if active
            x_minus, y_minus, z_minus = MOTION ? get_spin_coordinates(p_x, p_y, p_z, i, s_minus) : (x, y, z)
            x_plus, y_plus, z_plus = MOTION ? get_spin_coordinates(p_x, p_y, p_z, i, s_plus) : (x, y, z)

            Bx_minus, By_minus = reim(s_B1[s_minus])
            Bz_minus = x_minus * s_Gx[s_minus] + y_minus * s_Gy[s_minus] + z_minus * s_Gz[s_minus] + ΔBz - s_Δf[s_minus] * inv_γ
            Bx_plus, By_plus = reim(s_B1[s_plus])
            Bz_plus = x_plus * s_Gx[s_plus] + y_plus * s_Gy[s_plus] + z_plus * s_Gz[s_plus] + ΔBz - s_Δf[s_plus] * inv_γ

            Δt = s_Δt[s_idx] + s_Δt[s_minus] + s_Δt[s_plus]
            θx, θy, θz = rotation_vector(
                Bx_minus, By_minus, Bz_minus,
                Bx_plus, By_plus, Bz_plus,
                Δt,
                sim_method,
            )
            M_norm = mag_norm(T, Mxy_r, Mxy_i, Mz)
            Mxy_new_r, Mxy_new_i, Mz_new = rotate_magnetization(θx, θy, θz, Mxy_r, Mxy_i, Mz, T)
            Mxy_new_r, Mxy_new_i, Mz_new = restore_mag_norm(M_norm, Mxy_new_r, Mxy_new_i, Mz_new) # For reduced float precision only.

            E1 = exp(Δt * neg_inv_T1)
            E2 = exp(Δt * neg_inv_T2)
            Mxy_r = Mxy_new_r * E2
            Mxy_i = Mxy_new_i * E2
            Mz = Mz_new * E1 + ρ * (T(1) - E1)
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
