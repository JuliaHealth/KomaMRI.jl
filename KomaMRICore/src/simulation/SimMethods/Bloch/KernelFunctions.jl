using KernelAbstractions: @kernel, @Const, @index, @uniform, @groupsize, @localmem

## COV_EXCL_START

@kernel function apply_excitation!(Mxy, Mz, @Const(φ), @Const(B1), @Const(Bz), @Const(B), @Const(ΔT1), @Const(ΔT2), @Const(ρ))
    i_g = @index(Global)
    i_l = @index(Local)

    @uniform T = eltype(φ)
    @uniform N = @groupsize()[1]
    @uniform N_Δt = size(φ, 2)

    s_α_r = @localmem T (N,)
    s_α_i = @localmem T (N,)
    s_β_i = @localmem T (N,)
    s_β_r = @localmem T (N,)
    s_Mxy_r = @localmem T (N,)
    s_Mxy_i = @localmem T (N,)
    s_Mxy_new_r = @localmem T (N,)
    s_Mxy_new_i = @localmem T (N,)
    s_Mz = @localmem T (N,)
    s_Mz_new = @localmem T (N,)
    s_ρ = @localmem T (N,)

    @inbounds s_Mxy_r[i_l] = real(Mxy[i_g])
    @inbounds s_Mxy_i[i_l]  = imag(Mxy[i_g])
    @inbounds s_Mz[i_l] = Mz[i_g]
    @inbounds s_ρ[i_l] = ρ[i_g]

    @inbounds for t = 1 : N_Δt
        sin_φ, cos_φ = sincos(φ[i_g, t])
        s_α_r[i_l] = cos_φ
        if (iszero(B[i_g, t]))
            s_α_i[i_l] = -(Bz[i_g, t] / (B[i_g, t] + eps(T))) * sin_φ
            s_β_r[i_l] = (imag(B1[t]) / (B[i_g, t] + eps(T))) * sin_φ
            s_β_i[i_l] = -(real(B1[t]) / (B[i_g, t] + eps(T))) * sin_φ
        else
            s_α_i[i_l] = -(Bz[i_g, t] / B[i_g, t]) * sin_φ
            s_β_r[i_l] = (imag(B1[t]) / B[i_g, t]) * sin_φ
            s_β_i[i_l] = -(real(B1[t]) / B[i_g, t]) * sin_φ
        end
            s_Mxy_new_r[i_l] = 2 * (s_Mxy_i[i_l] * (s_α_r[i_l] * s_α_i[i_l] - s_β_r[i_l] * s_β_i[i_l]) +
                                s_Mz[i_l] * (s_α_i[i_l] * s_β_i[i_l] + s_α_r[i_l] * s_β_r[i_l])) +
                                s_Mxy_r[i_l] * (s_α_r[i_l]^2 - s_α_i[i_l]^2 - s_β_r[i_l]^2 + s_β_i[i_l]^2)
            s_Mxy_new_i[i_l] = -2 * (s_Mxy_r[i_l] * (s_α_r[i_l] * s_α_i[i_l] + s_β_r[i_l] * s_β_i[i_l]) -
                                s_Mz[i_l] * (s_α_r[i_l] * s_β_i[i_l] - s_α_i[i_l] * s_β_r[i_l])) +
                                s_Mxy_i[i_l] * (s_α_r[i_l]^2 - s_α_i[i_l]^2 + s_β_r[i_l]^2 - s_β_i[i_l]^2)
            s_Mz_new[i_l] = s_Mz[i_l] * (s_α_r[i_l]^2 + s_α_i[i_l]^2 - s_β_r[i_l]^2 - s_β_i[i_l]^2) -
                            2 * (s_Mxy_r[i_l] * (s_α_r[i_l] * s_β_r[i_l] - s_α_i[i_l] * s_β_i[i_l]) +
                            s_Mxy_i[i_l] * (s_α_r[i_l] * s_β_i[i_l] + s_α_i[i_l] * s_β_r[i_l]))
        s_Mxy_r[i_l] = s_Mxy_new_r[i_l] * ΔT2[i_g, t]
        s_Mxy_i[i_l] = s_Mxy_new_i[i_l] * ΔT2[i_g, t]
        s_Mz[i_l] = s_Mz_new[i_l] * ΔT1[i_g, t] + s_ρ[i_l] * (1 - ΔT1[i_g, t])
    end

    @inbounds Mxy[i_g] = s_Mxy_r[i_l] + 1im * s_Mxy_i[i_l]
    @inbounds Mz[i_g] = s_Mz[i_l]
end

## COV_EXCL_STOP