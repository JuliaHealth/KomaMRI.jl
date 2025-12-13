## COV_EXCL_START

@inline function effective_rotation_vector(
    Bx_prev,
    By_prev,
    Bz_prev,
    Bx_next,
    By_next,
    Bz_next,
    Δt,
    sim_method::SM,
) where {SM <: Union{Bloch,BlochMagnus1}}
    θx = Bx_prev * Δt
    θy = By_prev * Δt
    θz = Bz_prev * Δt
    return θx, θy, θz
end

@inline function effective_rotation_vector(
    Bx_prev, By_prev, Bz_prev, Bx_next, By_next, Bz_next, Δt, sim_method::BlochMagnus2
)
    τ = Δt / 2
    θx = (Bx_prev + Bx_next) * τ
    θy = (By_prev + By_next) * τ
    θz = (Bz_prev + Bz_next) * τ
    return θx, θy, θz
end

@inline function effective_rotation_vector(
    Bx_prev, By_prev, Bz_prev, Bx_next, By_next, Bz_next, Δt, sim_method::BlochMagnus4
)
    θx, θy, θz = effective_rotation_vector(
        Bx_prev, By_prev, Bz_prev, Bx_next, By_next, Bz_next, Δt, BlochMagnus2()
    )
    τ = Δt^2 / 12
    θx -= (By_prev * Bz_next - Bz_prev * By_next) * τ
    θy -= (Bz_prev * Bx_next - Bx_prev * Bz_next) * τ
    θz -= (Bx_prev * By_next - By_prev * Bx_next) * τ
    return θx, θy, θz
end

## COV_EXCL_STOP