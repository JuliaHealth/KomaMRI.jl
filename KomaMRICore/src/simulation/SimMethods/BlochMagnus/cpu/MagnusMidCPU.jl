function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::BlochMagnusMid2,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusMidCPUPrealloc
)
    T = eltype(p.ρ)
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_m, ωz_m, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    sample = 1

    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end

    @trace track_numbers=false for i in firstindex(seq.Δt):2:(lastindex(seq.Δt)-1)
        i_mid = i + 1
        i1 = i_mid + 1
        Δt = seq.t[i1] - seq.t[i]
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i_mid],
        )
        @. ωxy_m = seq.B1[i_mid] * B_to_ω
        @. ωz_m  = (seq.Gx[i_mid] * x + seq.Gy[i_mid] * y + seq.Gz[i_mid] * z + ΔBz) * B_to_ω + seq.Δf[i_mid] * T(2π)

        rotation_vector!(θxy, θz, ωxy_m, ωz_m, Δt, sim_method)
        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.

        @. M.xy = M.xy * exp(-Δt / p.T2)
        @. M.z = M.z * exp(-Δt / p.T1) + p.ρ * (T(1) - exp(-Δt / p.T1))
        outflow_spin_reset_at!(M, seq.t, i1, p.motion; replace_by=p.ρ)
        if !isempty(sig) && seq.ADC[i1]
            coords = spin_coordinates!(
                prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i1],
            )
            update_sensitivities!(prealloc.sens, sys.receiver, coords, p.motion)
            acquire_signal!(@view(sig[sample, :]), M.xy, prealloc.sens, coords)
            sample += 1
        end
    end

    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
