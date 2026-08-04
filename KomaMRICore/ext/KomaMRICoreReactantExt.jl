module KomaMRICoreReactantExt

using Reactant

import KomaMRIBase: DiscreteSequence, Phantom, γ
import KomaMRICore
import KomaMRICore: Bloch, BlochCPUPrealloc, Mag, PreallocResult, Spinor, mul!,
    outflow_spin_reset!, outflow_spin_reset_at!, spin_coordinates

const TracedArray = AbstractArray{<:Reactant.TracedRNumber}

function KomaMRICore.run_spin_precession!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::TracedArray,
    M::Mag,
    sim_method::Bloch,
    groupsize,
    backend::KomaMRICore.KA.CPU,
    prealloc::PreallocResult,
)
    T = eltype(p.ρ)
    Bz_old = prealloc.Bz_old
    Bz_new = prealloc.Bz_new
    ϕ = prealloc.ϕ
    Mxy = prealloc.M.xy
    ΔBz = prealloc.ΔBz
    fill!(ϕ, zero(T))
    block_time = zero(T)
    sample = 1
    x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[1])
    @. Bz_old = x * seq.Gx[1] + y * seq.Gy[1] + z * seq.Gz[1] + ΔBz
    for i in eachindex(seq.Δt)
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        @. Bz_new = x * seq.Gx[i + 1] + y * seq.Gy[i + 1] + z * seq.Gz[i + 1] + ΔBz
        @. ϕ += (Bz_old + Bz_new) * T(-π * γ) * seq.Δt[i]
        block_time += seq.Δt[i]
        if seq.ADC[i + 1]
            @. Mxy = exp(-block_time / p.T2) * M.xy * cis(ϕ)
            outflow_spin_reset!(Mxy, seq.t[i + 1], p.motion)
            Reactant.@allowscalar sig[sample] = sum(Mxy)
            sample += 1
        end
        Bz_old .= Bz_new
    end
    @. M.xy = M.xy * exp(-block_time / p.T2) * cis(ϕ)
    @. M.z = M.z * exp(-block_time / p.T1) + p.ρ * (T(1) - exp(-block_time / p.T1))
    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function KomaMRICore.run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::TracedArray,
    M::Mag,
    sim_method::Bloch,
    groupsize,
    backend::KomaMRICore.KA.CPU,
    prealloc::BlochCPUPrealloc,
)
    T = eltype(p.ρ)
    Bz = prealloc.Bz_old
    B = prealloc.Bz_new
    φ_half = prealloc.ϕ
    α = prealloc.Rot.α
    β = prealloc.Rot.β
    ΔBz = prealloc.ΔBz
    Maux_xy = prealloc.M.xy
    Maux_z = prealloc.M.z
    sample = 1
    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    for i in eachindex(seq.Δt)
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i])
        B1 = Reactant.@allowscalar seq.B1[i]
        @. Bz = (seq.Gx[i] * x + seq.Gy[i] * y + seq.Gz[i] * z) + ΔBz - seq.Δf[i] / T(γ)
        @. B = sqrt(abs2(B1) + Bz^2)
        @. φ_half = T(-π * γ) * (B * seq.Δt[i])
        @. α = cos(φ_half)
        @. B = T(-π * γ) * seq.Δt[i] * sinc(φ_half / T(π))
        @. α -= complex(zero(Bz), Bz * B)
        @. β = complex(imag(B1) * B, -real(B1) * B)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        outflow_spin_reset_at!(M, seq.t, i + 1, p.motion; replace_by=p.ρ)
        if seq.ADC[i + 1]
            Reactant.@allowscalar sig[sample] = sum(M.xy)
            sample += 1
        end
    end
    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end

end
