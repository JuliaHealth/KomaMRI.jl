# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_minus, ωz_minus, ωxy_plus, ωz_plus, Δt, sim_method::BlochMagnusGL2)
    @. θxy = (ωxy_minus + ωxy_plus) * (Δt / 2)
    @. θz  = (ωz_minus  + ωz_plus)  * (Δt / 2)
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_minus, By_minus, Bz_minus,
    Bx_plus, By_plus, Bz_plus,
    Δt,
    sim_method::BlochMagnusGL2,
)
    T = typeof(Δt)
    B_to_ω = T(-2π * γ)
    θ1_scale = B_to_ω * Δt / 2
    θx = (Bx_minus + Bx_plus) * θ1_scale
    θy = (By_minus + By_plus) * θ1_scale
    θz = (Bz_minus + Bz_plus) * θ1_scale
    return θx, θy, θz
end
