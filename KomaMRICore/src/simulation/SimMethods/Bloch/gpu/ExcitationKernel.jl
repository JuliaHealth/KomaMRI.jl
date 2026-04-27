## COV_EXCL_START

@kernel unsafe_indices=true inbounds=true function excitation_kernel!(
    sig_output::AbstractMatrix{Complex{T}}, 
    M_xy::AbstractVector{Complex{T}}, M_z, 
    @Const(p_x), @Const(p_y), @Const(p_z), @Const(p_ΔBz), @Const(p_T1), @Const(p_T2), @Const(p_ρ), N_spins,
    @Const(s_Gx), @Const(s_Gy), @Const(s_Gz), @Const(s_Δt), @Const(s_Δf), @Const(s_B1), @Const(s_ADC), s_length,
    ::Val{MOTION}, ::Val{USE_WARP_REDUCTION},
    sim_method::SM
) where {T, MOTION, USE_WARP_REDUCTION, SM <: BlochLikeSimMethods}

    @uniform N = @groupsize()[1]
    i_l = @index(Local, Linear)
    i_g = @index(Group, Linear)
    i = (i_g - 1u32) * UInt32(N) + i_l

    sig_group_r = @localmem T USE_WARP_REDUCTION ? 32 : N
    sig_group_i = @localmem T USE_WARP_REDUCTION ? 32 : N
    sig_r = zero(T)
    sig_i = zero(T)

    if i <= N_spins
        ΔBz = p_ΔBz[i]
        Mxy_r, Mxy_i = reim(M_xy[i])
        Mz = M_z[i]
        ρ = p_ρ[i]
        T1 = p_T1[i]
        T2 = p_T2[i]

        # Calculate initial B
        x, y, z = get_spin_coordinates(p_x, p_y, p_z, i, 1)
        Bx_prev, By_prev = reim(s_B1[1])
        Bz_prev = x * s_Gx[1] + y * s_Gy[1] + z * s_Gz[1] + ΔBz - s_Δf[1] / T(γ)

        ADC_idx = 1u32
        s_idx = 2u32
        while (s_idx <= s_length)
            if MOTION
                x, y, z = get_spin_coordinates(p_x, p_y, p_z, i, s_idx)
            end
            Bx_next, By_next = reim(s_B1[s_idx])
            Bz_next = (x * s_Gx[s_idx] + y * s_Gy[s_idx] + z * s_Gz[s_idx]) + ΔBz - s_Δf[s_idx] / T(γ)
            
            Δt = s_Δt[s_idx - 1]

            # Spinor rotation
            θx, θy, θz = effective_rotation_vector(Bx_prev, By_prev, Bz_prev, Bx_next, By_next, Bz_next, Δt, sim_method)
            θ =  sqrt(θx^2 + θy^2 + θz^2)
            sin_θ, cos_θ = sincos(T(-π * γ) * θ)
            α_r = cos_θ
            if iszero(θ)
                α_i = -(θz / (θ + eps(T))) * sin_θ
                β_r =  (θy / (θ + eps(T))) * sin_θ
                β_i = -(θx / (θ + eps(T))) * sin_θ
            else
                α_i = -(θz / θ) * sin_θ
                β_r =  (θy / θ) * sin_θ
                β_i = -(θx / θ) * sin_θ
            end

            Mxy_new_r = 2 * (Mxy_i * (α_r * α_i - β_r * β_i) +
                        Mz * (α_i * β_i + α_r * β_r)) +
                        Mxy_r * (α_r^2 - α_i^2 - β_r^2 + β_i^2)
            
            Mxy_new_i = -2 * (Mxy_r * (α_r * α_i + β_r * β_i) -
                        Mz * (α_r * β_i - α_i * β_r)) +
                        Mxy_i * (α_r^2 - α_i^2 + β_r^2 - β_i^2)
            
            Mz_new =    Mz * (α_r^2 + α_i^2 - β_r^2 - β_i^2) -
                        2 * (Mxy_r * (α_r * β_r - α_i * β_i) +
                        Mxy_i * (α_r * β_i + α_i * β_r))
            
            # Relaxation
            E1 = exp(-Δt / T1)
            E2 = exp(-Δt / T2)
            Mxy_r = Mxy_new_r * E2
            Mxy_i = Mxy_new_i * E2
            Mz = Mz_new * E1 + ρ * (T(1) - E1)

            # Acquire Signal
            if s_idx <= s_length && s_ADC[s_idx]
                sig_r, sig_i = reduce_signal!(Mxy_r, Mxy_i, sig_group_r, sig_group_i, i_l, N, T, Val(USE_WARP_REDUCTION))
                if i_l == 1u32
                    sig_output[i_g, ADC_idx] = complex(sig_r, sig_i)
                end
                ADC_idx += 1u32
            end

            Bx_prev, By_prev, Bz_prev = Bx_next, By_next, Bz_next
            s_idx += 1u32
        end

        M_xy[i] = complex(Mxy_r, Mxy_i)
        M_z[i] = Mz
    end
end

## COV_EXCL_STOP