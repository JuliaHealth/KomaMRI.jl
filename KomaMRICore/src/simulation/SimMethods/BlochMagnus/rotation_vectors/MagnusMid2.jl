# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_m, ωz_m, Δt, sim_method::BlochMagnusMid2)
    @. θxy = ωxy_m * Δt
    @. θz  = ωz_m  * Δt
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(
    Bx_m, By_m, Bz_m,
    Δt,
    B_to_ω,
    sim_method::BlochMagnusMid2,
)
    θ_scale = B_to_ω * Δt
    θx = Bx_m * θ_scale
    θy = By_m * θ_scale
    θz = Bz_m * θ_scale
    return θx, θy, θz
end

@inline function rotation_vector(
    Bx_m, By_m, Bz_m,
    Δt,
    sim_method::BlochMagnusMid2,
)
    B_to_ω = typeof(Δt)(-2π * γ)
    return rotation_vector(Bx_m, By_m, Bz_m, Δt, B_to_ω, sim_method)
end
