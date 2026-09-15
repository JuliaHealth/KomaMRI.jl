function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::Union{BlochMagnusQuad2,BlochMagnusQuad4},
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusQuadCPUPrealloc{T}
) where {T<:Real}
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_0, ωz_0, ωxy_m, ωz_m, ωxy_1, ωz_1, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    sample = 1

    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end

    i = firstindex(seq.Δt)
    cached_i0 = 0
    while i + 1 <= lastindex(seq.Δt)
        i0 = i
        im = i0 + 1
        i1 = im + 1
        (seq.excitation_bool[i0] && seq.excitation_bool[im]) ||
            throw(ArgumentError("BlochMagnusQuad RF intervals must contain midpoint samples."))

        Δt = seq.t[i1] - seq.t[i0]

        if cached_i0 != i0
            x, y, z = spin_coordinates!(
                prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i0],
            )
            @. ωxy_0 = seq.B1[i0] * B_to_ω
            @. ωz_0  = (seq.Gx[i0] * x + seq.Gy[i0] * y + seq.Gz[i0] * z + ΔBz) * B_to_ω + seq.Δf[i0] * T(2π)
        end
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[im],
        )
        @. ωxy_m = seq.B1[im] * B_to_ω
        @. ωz_m  = (seq.Gx[im] * x + seq.Gy[im] * y + seq.Gz[im] * z + ΔBz) * B_to_ω + seq.Δf[im] * T(2π)
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i1],
        )
        @. ωxy_1 = seq.B1[i1] * B_to_ω
        @. ωz_1  = (seq.Gx[i1] * x + seq.Gy[i1] * y + seq.Gz[i1] * z + ΔBz) * B_to_ω + seq.Δf[i1] * T(2π)

        rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_m, ωz_m, ωxy_1, ωz_1, Δt, sim_method)

        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.

        @. M.xy = M.xy * exp(-Δt / p.T2)
        @. M.z = M.z * exp(-Δt / p.T1) + p.ρ * (T(1) - exp(-Δt / p.T1))
        outflow_spin_reset_at!(M, seq.t, i1, p.motion; replace_by=p.ρ)
        if seq.ADC[i1]
            update_sensitivities!(
                prealloc.sens, sys.receiver, (x, y, z), p.motion,
            )
            acquire_signal!(
                @view(sig[sample, :]), M.xy, prealloc.sens, (x, y, z),
            )
            sample += 1
        end
        cached_i0 = i1
        ωxy_0, ωxy_1 = ωxy_1, ωxy_0
        ωz_0, ωz_1 = ωz_1, ωz_0
        i = i1
    end

    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
