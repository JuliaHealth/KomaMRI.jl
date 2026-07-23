# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_m, ωz_m, ωxy_1, ωz_1, Δt, sim_method::BlochMagnusQuad2)
    @. θxy = (ωxy_0 + 4ωxy_m + ωxy_1) * (Δt / 6)
    @. θz  = (ωz_0  + 4ωz_m  + ωz_1)  * (Δt / 6)
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_0, By_0, Bz_0,
    Bx_m, By_m, Bz_m,
    Bx_1, By_1, Bz_1,
    Δt,
    sim_method::BlochMagnusQuad2,
)
    B_to_ω = typeof(Δt)(-2π * γ)
    θ1_scale = B_to_ω * Δt / 6
    θx = (Bx_0 + 4Bx_m + Bx_1) * θ1_scale
    θy = (By_0 + 4By_m + By_1) * θ1_scale
    θz = (Bz_0 + 4Bz_m + Bz_1) * θ1_scale
    return θx, θy, θz
end
