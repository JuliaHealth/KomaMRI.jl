# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_0, ωz_0, Δt, sim_method::BlochMagnusConst1)
    @. θxy = ωxy_0 * Δt
    @. θz  = ωz_0 * Δt
    return nothing
end

function rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_1, ωz_1, Δt, sim_method::BlochMagnusConst1)
    return rotation_vector!(θxy, θz, ωxy_0, ωz_0, Δt, sim_method)
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(Bx_0, By_0, Bz_0, Δt, sim_method::SM) where {SM<:Union{Bloch,BlochMagnusConst1}}
    B_to_ω = typeof(Δt)(-2π * γ)
    θx = Bx_0 * B_to_ω * Δt
    θy = By_0 * B_to_ω * Δt
    θz = Bz_0 * B_to_ω * Δt
    return θx, θy, θz
end

@inline function rotation_vector(Bx_0, By_0, Bz_0, Bx_1, By_1, Bz_1, Δt, sim_method::SM) where {SM<:Union{Bloch,BlochMagnusConst1}}
    return rotation_vector(Bx_0, By_0, Bz_0, Δt, sim_method)
end
