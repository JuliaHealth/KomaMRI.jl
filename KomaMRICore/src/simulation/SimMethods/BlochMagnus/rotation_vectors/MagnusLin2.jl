# -- CPU vectorized rotation vector ------------------------------------------
function rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_1, ωz_1, Δt, sim_method::BlochMagnusLin2)
    @. θxy = (ωxy_0 + ωxy_1) * (Δt / 2)
    @. θz  = (ωz_0 + ωz_1) * (Δt / 2)
    return nothing
end

# -- GPU scalar rotation vector ----------------------------------------------
@inline function rotation_vector(Bx_0, By_0, Bz_0, Bx_1, By_1, Bz_1, Δt, sim_method::BlochMagnusLin2)
    B_to_ω = typeof(Δt)(-2π * γ)
    half_Δt = Δt / 2
    θx = (Bx_0 + Bx_1) * B_to_ω * half_Δt
    θy = (By_0 + By_1) * B_to_ω * half_Δt
    θz = (Bz_0 + Bz_1) * B_to_ω * half_Δt
    return θx, θy, θz
end
