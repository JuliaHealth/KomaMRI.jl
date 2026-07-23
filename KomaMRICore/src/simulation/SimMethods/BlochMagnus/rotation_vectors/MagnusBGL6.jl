# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(
    θxy, θz,
    ωxy_minus, ωz_minus,
    ωxy_0, ωz_0,
    ωxy_plus, ωz_plus,
    i0xy, i0z,
    i1xy, i1z,
    i2xy, i2z,
    jxy, jz,
    boxxy, boxz,
    Δt,
    sim_method::BlochMagnusBGL6,
)
    T = typeof(Δt)
    sqrt15 = sqrt(T(15))
    Δt2 = Δt^2
    θ34_scale = T(3) / 5 * Δt2 * Δt

    @. i0xy = (5ωxy_minus + 8ωxy_0 + 5ωxy_plus) / 18
    @. i0z = (5ωz_minus + 8ωz_0 + 5ωz_plus) / 18
    @. i1xy = sqrt15 * (ωxy_plus - ωxy_minus) / 36
    @. i1z = sqrt15 * (ωz_plus - ωz_minus) / 36
    @. i2xy = (ωxy_minus + ωxy_plus) / 24
    @. i2z = (ωz_minus + ωz_plus) / 24

    @. jxy = T(3) / 2 * i0xy - T(6) * i2xy
    @. jz = T(3) / 2 * i0z - T(6) * i2z
    @. θxy = -im * Δt2 * (i1xy * jz - jxy * i1z)
    @. θz = Δt2 * imag(conj(i1xy) * jxy)

    @. boxxy = Δt / 2 * i2xy - θxy / 60
    @. boxz = Δt / 2 * i2z - θz / 60
    @. θxy = Δt * i0xy + θxy +
             Δt2 * (
                 i0xy * (real(conj(i0xy) * boxxy) + i0z * boxz) -
                 boxxy * (abs2(i0xy) + i0z^2)
             ) +
             θ34_scale * (
                 i1xy * (real(conj(i1xy) * jxy) + i1z * jz) -
                 jxy * (abs2(i1xy) + i1z^2)
             )
    @. θz = Δt * i0z + θz +
            Δt2 * (
                i0z * (real(conj(i0xy) * boxxy) + i0z * boxz) -
                boxz * (abs2(i0xy) + i0z^2)
            ) +
            θ34_scale * (
                i1z * (real(conj(i1xy) * jxy) + i1z * jz) -
                jz * (abs2(i1xy) + i1z^2)
            )
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_minus, By_minus, Bz_minus,
    Bx_0, By_0, Bz_0,
    Bx_plus, By_plus, Bz_plus,
    Δt,
    sim_method::BlochMagnusBGL6,
)
    T = typeof(Δt)
    B_to_ω = T(-2π * γ)
    ωx_minus, ωy_minus, ωz_minus = Bx_minus * B_to_ω, By_minus * B_to_ω, Bz_minus * B_to_ω
    ωx_0, ωy_0, ωz_0 = Bx_0 * B_to_ω, By_0 * B_to_ω, Bz_0 * B_to_ω
    ωx_plus, ωy_plus, ωz_plus = Bx_plus * B_to_ω, By_plus * B_to_ω, Bz_plus * B_to_ω

    i0x = (5ωx_minus + 8ωx_0 + 5ωx_plus) / 18
    i0y = (5ωy_minus + 8ωy_0 + 5ωy_plus) / 18
    i0z = (5ωz_minus + 8ωz_0 + 5ωz_plus) / 18
    i1x = sqrt(T(15)) * (ωx_plus - ωx_minus) / 36
    i1y = sqrt(T(15)) * (ωy_plus - ωy_minus) / 36
    i1z = sqrt(T(15)) * (ωz_plus - ωz_minus) / 36
    i2x = (ωx_minus + ωx_plus) / 24
    i2y = (ωy_minus + ωy_plus) / 24
    i2z = (ωz_minus + ωz_plus) / 24

    jx = T(3) / 2 * i0x - T(6) * i2x
    jy = T(3) / 2 * i0y - T(6) * i2y
    jz = T(3) / 2 * i0z - T(6) * i2z
    Δt2 = Δt^2
    θ2x = Δt2 * (i1y * jz - i1z * jy)
    θ2y = Δt2 * (i1z * jx - i1x * jz)
    θ2z = Δt2 * (i1x * jy - i1y * jx)

    boxx = Δt / 2 * i2x - θ2x / 60
    boxy = Δt / 2 * i2y - θ2y / 60
    boxz = Δt / 2 * i2z - θ2z / 60
    dot0 = i0x * boxx + i0y * boxy + i0z * boxz
    norm0 = i0x^2 + i0y^2 + i0z^2
    dot1 = i1x * jx + i1y * jy + i1z * jz
    norm1 = i1x^2 + i1y^2 + i1z^2

    θ34_scale = T(3) / 5 * Δt2 * Δt
    θx = Δt * i0x + θ2x +
         Δt2 * (i0x * dot0 - boxx * norm0) +
         θ34_scale * (i1x * dot1 - jx * norm1)
    θy = Δt * i0y + θ2y +
         Δt2 * (i0y * dot0 - boxy * norm0) +
         θ34_scale * (i1y * dot1 - jy * norm1)
    θz = Δt * i0z + θ2z +
         Δt2 * (i0z * dot0 - boxz * norm0) +
         θ34_scale * (i1z * dot1 - jz * norm1)
    return θx, θy, θz
end
