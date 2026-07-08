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
    B_to_ω,
    B_to_ω2_sqrt3,
    sim_method::BlochMagnusGL2,
)
    θ1_scale = B_to_ω * Δt / 2
    θx = (Bx_minus + Bx_plus) * θ1_scale
    θy = (By_minus + By_plus) * θ1_scale
    θz = (Bz_minus + Bz_plus) * θ1_scale
    return θx, θy, θz
end

@inline function rotation_vector(
    Bx_minus, By_minus, Bz_minus,
    Bx_plus, By_plus, Bz_plus,
    Δt,
    sim_method::BlochMagnusGL2,
)
    B_to_ω = typeof(Δt)(-2π * γ)
    return rotation_vector(
        Bx_minus, By_minus, Bz_minus,
        Bx_plus, By_plus, Bz_plus,
        Δt,
        B_to_ω,
        B_to_ω * B_to_ω * sqrt(typeof(Δt)(3)),
        sim_method,
    )
end
