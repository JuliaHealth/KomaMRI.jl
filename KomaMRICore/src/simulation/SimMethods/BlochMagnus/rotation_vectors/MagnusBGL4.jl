# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(
    θxy, θz,
    ωxy_minus, ωz_minus,
    ωxy_center, ωz_center,
    ωxy_plus, ωz_plus,
    i0xy, i0z,
    i1xy, i1z,
    i2xy, i2z,
    jxy, jz,
    boxxy, boxz,
    Δt,
    sim_method::BlochMagnusBGL4,
)
    @. i0xy = (5ωxy_minus + 8ωxy_center + 5ωxy_plus) / 18
    @. i0z = (5ωz_minus + 8ωz_center + 5ωz_plus) / 18
    @. i1xy = sqrt(typeof(Δt)(15)) * (ωxy_plus - ωxy_minus) / 36
    @. i1z = sqrt(typeof(Δt)(15)) * (ωz_plus - ωz_minus) / 36

    Δt2 = Δt^2
    @. θxy = Δt * i0xy - im * Δt2 * (i1xy * i0z - i0xy * i1z)
    @. θz = Δt * i0z + Δt2 * imag(conj(i1xy) * i0xy)
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_minus, By_minus, Bz_minus,
    Bx_center, By_center, Bz_center,
    Bx_plus, By_plus, Bz_plus,
    Δt,
    sim_method::BlochMagnusBGL4,
)
    T = typeof(Δt)
    B_to_ω = T(-2π * γ)
    θ1_scale = B_to_ω * Δt / 18
    θ2_scale = B_to_ω^2 * sqrt(T(15)) * Δt^2 / 648

    Sx = 5Bx_minus + 8Bx_center + 5Bx_plus
    Sy = 5By_minus + 8By_center + 5By_plus
    Sz = 5Bz_minus + 8Bz_center + 5Bz_plus
    dBx = Bx_plus - Bx_minus
    dBy = By_plus - By_minus
    dBz = Bz_plus - Bz_minus

    θx = Sx * θ1_scale
    θy = Sy * θ1_scale
    θz = Sz * θ1_scale
    θx -= θ2_scale * (Sy * dBz - Sz * dBy)
    θy -= θ2_scale * (Sz * dBx - Sx * dBz)
    θz -= θ2_scale * (Sx * dBy - Sy * dBx)
    return θx, θy, θz
end
