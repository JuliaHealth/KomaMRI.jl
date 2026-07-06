# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_minus, ωz_minus, ωxy_plus, ωz_plus, Δt, sim_method::BlochMagnusGL4)
    @. θxy = (ωxy_minus + ωxy_plus) * (Δt / 2)
    @. θz  = (ωz_minus  + ωz_plus)  * (Δt / 2)
    θ2_scale = sqrt(typeof(Δt)(3)) * Δt^2 / 12
    @. θxy -= im * (ωxy_plus * ωz_minus - ωxy_minus * ωz_plus) * θ2_scale
    @. θz  += imag(conj(ωxy_plus) * ωxy_minus) * θ2_scale
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_minus, By_minus, Bz_minus,
    Bx_plus, By_plus, Bz_plus,
    Δt,
    sim_method::BlochMagnusGL4,
)
    T = typeof(Δt)
    B_to_ω = T(-2π * γ)
    θ1_scale = B_to_ω * Δt / 2
    θ2_scale = B_to_ω^2 * sqrt(T(3)) * Δt^2 / 12
    θx = (Bx_minus + Bx_plus) * θ1_scale
    θy = (By_minus + By_plus) * θ1_scale
    θz = (Bz_minus + Bz_plus) * θ1_scale
    θx += θ2_scale * (By_plus * Bz_minus - Bz_plus * By_minus)
    θy += θ2_scale * (Bz_plus * Bx_minus - Bx_plus * Bz_minus)
    θz += θ2_scale * (Bx_plus * By_minus - By_plus * Bx_minus)
    return θx, θy, θz
end
