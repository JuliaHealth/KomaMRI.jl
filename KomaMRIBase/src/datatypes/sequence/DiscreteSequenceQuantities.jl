# ===========================================================================
# 7. Sampled Sequence Derived Quantities
# ===========================================================================
#
# Quantities computed from an already-sampled sequence: k-space, gradient
# moments, slew rate, and eddy-current estimates.

# -- 7.1. Event waveform accessors ------------------------------------------
get_grads(seqd::DiscreteSequence) = seqd.Gx, seqd.Gy, seqd.Gz

get_rfs(seqd::DiscreteSequence) = seqd.B1, seqd.Δf, seqd.ψ

# -- 7.2. ADC-time sampling --------------------------------------------------
function values_at_adc_times(seqd, values, t)
    t_adc = seqd.t[seqd.ADC]
    out = zeros(eltype(values), length(t_adc), size(values, 2))
    for i in axes(values, 2)
        out[:, i] = linear_interpolate_samples((t=t, A=view(values, :, i)), t_adc; default=zero(eltype(values)))
    end
    return out
end

# -- 7.3. Gradient moments and k-space --------------------------------------
function get_Mk(seqd::DiscreteSequence, k; rf_idx=Int[], rf_types=RFUse[])
    get_sign(::Excitation) =  0
    get_sign(::Refocusing) = -1
    get_sign(::RFUse)      =  1
    t, Δt = seqd.t[1:end-1], seqd.Δt
    G = (seqd.Gx, seqd.Gy, seqd.Gz)
    mk = zeros(length(t), 3)
    parts = kfoldperm(length(t), 1; breaks=rf_idx)
    for i = 1:3
        mkf = 0
        for (rf, p) in enumerate(parts)
            mk[p,i] = cumtrapz(Δt[p]', [t[p]' t[p[end]]'.+Δt[p[end]]].^k .* G[i][p[1]:p[end]+1]')[:]
            rf > 1 && (mk[p,i] .+= mkf * get_sign(rf_types[rf-1]))
            mkf = mk[p[end],i]
        end
    end
    Mk = γ * mk
    return Mk, values_at_adc_times(seqd, Mk, seqd.t[2:end])
end
get_kspace(seqd::DiscreteSequence; kwargs...) = get_Mk(seqd, 0; kwargs...)
get_M1(seqd::DiscreteSequence; kwargs...) = get_Mk(seqd, 1; kwargs...)
get_M2(seqd::DiscreteSequence; kwargs...) = get_Mk(seqd, 2; kwargs...)

# -- 7.4. Slew rate ----------------------------------------------------------
function get_slew_rate(seqd::DiscreteSequence)
    Gx, Gy, Gz = get_grads(seqd)
    SR = [Gx[2:end] .- Gx[1:end-1] Gy[2:end] .- Gy[1:end-1] Gz[2:end] .- Gz[1:end-1]] ./ seqd.Δt
    SR[isnan.(SR)] .= 0.0
    length(seqd.Δt) >= 1 && (SR[1, :] .= 0.0)
    length(seqd.Δt) >= 2 && (SR[end, :] .= 0.0)
    return SR, values_at_adc_times(seqd, SR, seqd.t[2:end])
end

# -- 7.5. Eddy currents ------------------------------------------------------
function get_eddy_currents(seqd::DiscreteSequence; λ=80e-3)
    Gx, Gy, Gz = get_grads(seqd)
    t = seqd.t[1:end-1]
    dG = [Gx[2:end] .- Gx[1:end-1] Gy[2:end] .- Gy[1:end-1] Gz[2:end] .- Gz[1:end-1]]
    ec(t, λ) = exp.(-t ./ λ) .* (t .>= 0)
    EC = [sum(dG[:, j] .* ec(t[i] .- t, λ)) for i in eachindex(t), j = 1:3]
    return EC, values_at_adc_times(seqd, EC, seqd.t[2:end])
end
