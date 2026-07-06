function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sim_method::Union{BlochMagnusBGL4,BlochMagnusBGL6},
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusBGLCPUPrealloc
)
    T = eltype(p.ρ)
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_minus, ωz_minus, ωxy_center, ωz_center, ωxy_plus, ωz_plus,
        i0xy, i0z, i1xy, i1z, i2xy, i2z, jxy, jz, boxxy, boxz,
        θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    sample = 1

    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end

    i = firstindex(seq.Δt)
    while i + 3 <= lastindex(seq.Δt)
        i0 = i
        i_minus = i0 + 1
        i_center = i_minus + 1
        i_plus = i_center + 1
        i1 = i_plus + 1
        (seq.excitation_bool[i0] && seq.excitation_bool[i_minus] && seq.excitation_bool[i_center] && seq.excitation_bool[i_plus]) ||
            throw(ArgumentError("$(typeof(sim_method)) RF intervals must contain three Blanes Gauss-Legendre nodes."))

        Δt = seq.t[i1] - seq.t[i0]
        x_minus, y_minus, z_minus = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i_minus])
        x_center, y_center, z_center = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i_center])
        x_plus, y_plus, z_plus = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i_plus])

        @. ωxy_minus = seq.B1[i_minus] * B_to_ω
        @. ωz_minus = (seq.Gx[i_minus] * x_minus + seq.Gy[i_minus] * y_minus + seq.Gz[i_minus] * z_minus + ΔBz) * B_to_ω + seq.Δf[i_minus] * T(2π)
        @. ωxy_center = seq.B1[i_center] * B_to_ω
        @. ωz_center = (seq.Gx[i_center] * x_center + seq.Gy[i_center] * y_center + seq.Gz[i_center] * z_center + ΔBz) * B_to_ω + seq.Δf[i_center] * T(2π)
        @. ωxy_plus = seq.B1[i_plus] * B_to_ω
        @. ωz_plus = (seq.Gx[i_plus] * x_plus + seq.Gy[i_plus] * y_plus + seq.Gz[i_plus] * z_plus + ΔBz) * B_to_ω + seq.Δf[i_plus] * T(2π)

        rotation_vector!(
            θxy, θz,
            ωxy_minus, ωz_minus,
            ωxy_center, ωz_center,
            ωxy_plus, ωz_plus,
            i0xy, i0z, i1xy, i1z, i2xy, i2z, jxy, jz, boxxy, boxz,
            Δt,
            sim_method,
        )

        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.

        @. M.xy = M.xy * exp(-Δt / p.T2)
        @. M.z = M.z * exp(-Δt / p.T1) + p.ρ * (T(1) - exp(-Δt / p.T1))
        outflow_spin_reset_at!(M, seq.t, i1, p.motion; replace_by=p.ρ)
        if seq.ADC[i1]
            sig[sample] = sum(M.xy)
            sample += 1
        end
        i = i1
    end

    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
