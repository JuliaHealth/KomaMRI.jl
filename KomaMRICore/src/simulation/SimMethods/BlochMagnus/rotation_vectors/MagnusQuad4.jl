# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_m, ωz_m, ωxy_1, ωz_1, Δt, sim_method::BlochMagnusQuad4)
    @. θxy = (ωxy_0 + 4ωxy_m + ωxy_1) * (Δt / 6)
    @. θz  = (ωz_0  + 4ωz_m  + ωz_1)  * (Δt / 6)
    Δt2 = Δt^2
    @. θxy -= im * Δt2 * (
        (ωxy_1 * ωz_0 - ωxy_0 * ωz_1) / 12 +
        (
            (ωxy_1 - ωxy_0) * (ωz_m - (ωz_0 + ωz_1) / 2) -
            (ωxy_m - (ωxy_0 + ωxy_1) / 2) * (ωz_1 - ωz_0)
        ) / 15
    )
    @. θz += Δt2 * (
        imag(conj(ωxy_1) * ωxy_0) / 12 +
        imag(conj(ωxy_1 - ωxy_0) * (ωxy_m - (ωxy_0 + ωxy_1) / 2)) / 15
    )
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_0, By_0, Bz_0,
    Bx_m, By_m, Bz_m,
    Bx_1, By_1, Bz_1,
    Δt,
    sim_method::BlochMagnusQuad4,
)
    B_to_ω = typeof(Δt)(-2π * γ)
    θ1_scale = B_to_ω * Δt / 6
    θ2_scale = B_to_ω^2 * Δt^2
    θx = (Bx_0 + 4Bx_m + Bx_1) * θ1_scale
    θy = (By_0 + 4By_m + By_1) * θ1_scale
    θz = (Bz_0 + 4Bz_m + Bz_1) * θ1_scale

    δBx = Bx_m - (Bx_0 + Bx_1) / 2
    δBy = By_m - (By_0 + By_1) / 2
    δBz = Bz_m - (Bz_0 + Bz_1) / 2
    dBx = Bx_1 - Bx_0
    dBy = By_1 - By_0
    dBz = Bz_1 - Bz_0

    θx += θ2_scale * ((By_1 * Bz_0 - Bz_1 * By_0) / 12 + (dBy * δBz - dBz * δBy) / 15)
    θy += θ2_scale * ((Bz_1 * Bx_0 - Bx_1 * Bz_0) / 12 + (dBz * δBx - dBx * δBz) / 15)
    θz += θ2_scale * ((Bx_1 * By_0 - By_1 * Bx_0) / 12 + (dBx * δBy - dBy * δBx) / 15)
    return θx, θy, θz
end
