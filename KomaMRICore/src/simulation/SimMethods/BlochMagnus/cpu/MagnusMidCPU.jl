function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochMagnusMid2,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusMidCPUPrealloc{T}
) where {T<:Real}
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_m, ωz_m, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    (; neg_inv_T1, neg_inv_T2, E1) = prealloc.relaxation
    sample = 1

    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end

    i = firstindex(seq.Δt)
    while i + 1 <= lastindex(seq.Δt)
        im = i + 1
        i1 = im + 1
        Δt = seq.t[i1] - seq.t[i]
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[im])
        @. ωxy_m = seq.B1[im] * B_to_ω
        @. ωz_m  = (seq.Gx[im] * x + seq.Gy[im] * y + seq.Gz[im] * z + ΔBz) * B_to_ω + seq.Δf[im] * T(2π)

        rotation_vector!(θxy, θz, ωxy_m, ωz_m, Δt, sim_method)
        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.

        @. M.xy = M.xy * exp(Δt * neg_inv_T2)
        @. E1 = exp(Δt * neg_inv_T1)
        @. M.z = M.z * E1 + p.ρ * (T(1) - E1)
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
