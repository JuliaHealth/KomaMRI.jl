function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::Union{BlochMagnusGL2,BlochMagnusGL4},
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusGLCPUPrealloc
)
    T = eltype(p.ρ)
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_minus, ωz_minus, ωxy_plus, ωz_plus, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    sample = _sample_counter(sig)

    ψ_start = seq.ψ[1]
    @trace track_numbers=false if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end

    @trace track_numbers=false for i in firstindex(seq.Δt):3:(lastindex(seq.Δt)-2)
        i0 = i
        i_minus = i0 + 1
        i_plus = i_minus + 1
        i1 = i_plus + 1

        Δt = seq.t[i1] - seq.t[i0]
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i_minus],
        )
        @. ωxy_minus = seq.B1[i_minus] * B_to_ω
        @. ωz_minus = (seq.Gx[i_minus] * x + seq.Gy[i_minus] * y + seq.Gz[i_minus] * z + ΔBz) * B_to_ω + seq.Δf[i_minus] * T(2π)
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i_plus],
        )
        @. ωxy_plus = seq.B1[i_plus] * B_to_ω
        @. ωz_plus = (seq.Gx[i_plus] * x + seq.Gy[i_plus] * y + seq.Gz[i_plus] * z + ΔBz) * B_to_ω + seq.Δf[i_plus] * T(2π)

        rotation_vector!(θxy, θz, ωxy_minus, ωz_minus, ωxy_plus, ωz_plus, Δt, sim_method)

        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.

        @. M.xy = M.xy * exp(-Δt / p.T2)
        @. M.z = M.z * exp(-Δt / p.T1) + p.ρ * (T(1) - exp(-Δt / p.T1))
        outflow_spin_reset_at!(M, seq.t, i1, p.motion; replace_by=p.ρ)
        if !isempty(sig)
            adc = seq.ADC[i1]
            acquire_signal!(sig, sample, adc, M, sys, prealloc, p, seq, i1 - 1)
            sample += adc
        end
    end

    ψ_end = seq.ψ[end]
    @trace track_numbers=false if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
