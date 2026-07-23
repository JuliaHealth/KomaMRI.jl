# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_1, ωz_1, Δt, sim_method::BlochMagnusLinComm2)
    @. θxy = (ωxy_0 + ωxy_1) * (Δt / 2)
    @. θz  = (ωz_0 + ωz_1) * (Δt / 2)
    θ2_scale = Δt^2 / 12
    @. θxy -= im * (ωxy_1 * ωz_0 - ωxy_0 * ωz_1) * θ2_scale
    @. θz  += imag(conj(ωxy_1) * ωxy_0) * θ2_scale
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(Bx_0, By_0, Bz_0, Bx_1, By_1, Bz_1, Δt, sim_method::BlochMagnusLinComm2)
    B_to_ω = typeof(Δt)(-2π * γ)
    θ1_scale = B_to_ω * Δt / 2
    θ2_scale = B_to_ω^2 * Δt^2 / 12
    θx = (Bx_0 + Bx_1) * θ1_scale
    θy = (By_0 + By_1) * θ1_scale
    θz = (Bz_0 + Bz_1) * θ1_scale
    θx += θ2_scale * (By_1 * Bz_0 - Bz_1 * By_0)
    θy += θ2_scale * (Bz_1 * Bx_0 - Bx_1 * Bz_0)
    θz += θ2_scale * (Bx_1 * By_0 - By_1 * Bx_0)
    return θx, θy, θz
end
