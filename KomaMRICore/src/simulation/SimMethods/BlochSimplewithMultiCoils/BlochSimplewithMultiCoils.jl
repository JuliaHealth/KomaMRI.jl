#Simplest sim method, works for GPU and CPU but not optimized for either. Although Bloch()
#is the simulation method chosen if none is passed, the run_spin_precession! and
#run_spin_excitation! functions in this file are dispatched to at the most abstract level,
#so new simulation methods will start by using these functions.
struct BlochSimplewithMultiCoils{T<:Real} <: SimulationMethod
    ncoils::Int
    radius::T
    L::T
end

BlochSimplewithMultiCoils(; ncoils=8, radius=0.20, L=1.00) = 
    BlochSimplewithMultiCoils(ncoils, radius, L)

export BlochSimplewithMultiCoils

function sim_output_dim(
    obj::Phantom, seq::Sequence, sys::Scanner, sim_method::BlochSimplewithMultiCoils
)
    return (sum(seq.ADC.N), sim_method.ncoils)
end

"""
    run_spin_precession(obj, seq, Xt, sig)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
precession.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochSimplewithMultiCoils,
    groupsize,
    backend::KA.Backend,
    prealloc::PreallocResult
) where {T<:Real}
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    #Effective field
    Bz = x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz' .+ p.Δw ./ T(2π .* γ)
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π .* γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π .* γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy precession and relaxation, and Mz relaxation
    tp   = cumsum(seq.Δt) # t' = t - t0
    dur  = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy  = M.xy .* exp.(-tp' ./ p.T2) .* cis.(ϕ) #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(Mxy, seq.t[2:end]', p.motion)
    outflow_spin_reset!(M, seq.t[2:end]', p.motion; replace_by=p.ρ)
    #Acquired signal
    acquire_signal_precession!(sig, Mxy, seq, p.motion, x, y, z, sim_method)
    return nothing
end

function birdcage_sensitivities(
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T};
    ncoils=8,
    radius=T(0.20),
    L=T(0.30)
) where {T<:Real}

    nspins = length(x)
    sens = Matrix{Complex{T}}(undef, ncoils, nspins)

    ϵ = eps(T)

    for n in 1:ncoils
        ϕn = T(2π) * T(n - 1) / T(ncoils)

        xn = radius * cos(ϕn)
        yn = radius * sin(ϕn)

        for j in 1:nspins
            rn = sqrt((x[j] - xn)^2 + (y[j] - yn)^2 + ϵ)

            Bϕ =
                (1 / rn) *
                (
                    (L - z[j]) / sqrt(rn^2 + (L - z[j])^2) +
                    (L + z[j]) / sqrt(rn^2 + (L + z[j])^2)
                )

            Bx = -Bϕ * (y[j] - yn) / rn
            By =  Bϕ * (x[j] - xn) / rn

            sens[n, j] = Complex{T}(0.5) * (Bx - im * By)
        end

        sens[n, :] ./= maximum(abs.(sens[n, :])) + ϵ
    end

    return sens
end

"""
    M0 = run_spin_excitation(obj, seq, M0)

It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (actually, it's
    a part of the complete Mag vector and it's a part of the initial state for the next
    precession simulation step)
"""
function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochSimplewithMultiCoils,
    groupsize,
    backend::KA.Backend,
    prealloc::PreallocResult
) where {T<:Real}
    sample = 1
    sens = static_sensitivities(p, p.motion, sim_method)
    # Rotating frame -> RF frame
    ψ_start = @view seq.ψ[1:1]
    @. M.xy = M.xy * cis(-ψ_start)
    #Simulation
    for i in eachindex(seq.Δt)
        s = @views ( # This was the previous behaviour of seq[i], but it was hidden
            t = seq.t[i, :], tnew = seq.t[i + 1, :], Δt = seq.Δt[i, :],
            Gx = seq.Gx[i, :], Gy = seq.Gy[i, :], Gz = seq.Gz[i, :],
            B1 = seq.B1[i, :], Δf = seq.Δf[i, :],
            ADC = any(seq.ADC[i + 1, :])
        )
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, s.t)
        #Effective field
        ΔBz = p.Δw ./ T(2π .* γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        Bz = (s.Gx .* x .+ s.Gy .* y .+ s.Gz .* z) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B .+= (B .== 0) .* eps(T)
        #Spinor Rotation
        φ = T(-2π .* γ) .* (B .* s.Δt)
        mul!(Q(φ, s.B1 ./ B, Bz ./ B), M)
        #Relaxation
        @. M.xy = M.xy * exp(-s.Δt / p.T2)
        @. M.z  = M.z * exp(-s.Δt / p.T1) + p.ρ * (1 - exp(-s.Δt / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset!(M, s.tnew, p.motion; replace_by=p.ρ)
        #Acquire signal
        if s.ADC # ADC at the end of the time step
            acquire_signal!(sig, sample, M, p.motion, x, y, z, sim_method, sens)
            sample += 1
        end
    end
    # RF frame -> Rotating frame
    ψ_end = @view seq.ψ[end:end]
    @. M.xy = M.xy * cis(ψ_end)
    return nothing
end

function static_sensitivities(
    p::Phantom{T}, ::NoMotion, sim_method::BlochSimplewithMultiCoils
) where {T<:Real}
    return birdcage_sensitivities(
        p.x, p.y, p.z;
        ncoils=sim_method.ncoils,
        radius=T(sim_method.radius),
        L=T(sim_method.L)
    )
end

static_sensitivities(p::Phantom, ::Union{Motion, MotionList}, sim_method::BlochSimplewithMultiCoils) = nothing

function acquire_signal_precession!(
    sig,
    Mxy,
    seq::DiscreteSequence{T},
    ::NoMotion,
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    sim_method::BlochSimplewithMultiCoils,
) where {T<:Real}
    sens = birdcage_sensitivities(
        x, y, z;
        ncoils=sim_method.ncoils,
        radius=T(sim_method.radius),
        L=T(sim_method.L)
    )
    sig .= transpose(sens * (@view Mxy[:, findall(seq.ADC[2:end])]))
end

function acquire_signal_precession!(
    sig,
    Mxy,
    seq::DiscreteSequence{T},
    ::Union{Motion{T}, MotionList{T}},
    x::AbstractMatrix{T},
    y::AbstractMatrix{T},
    z::AbstractMatrix{T},
    sim_method::BlochSimplewithMultiCoils,
) where {T<:Real}
    adc = findall(seq.ADC[2:end])
    for (sample, i) in pairs(adc)
        sens = birdcage_sensitivities(
            @view(x[:, i]), @view(y[:, i]), @view(z[:, i]);
            ncoils=sim_method.ncoils,
            radius=T(sim_method.radius),
            L=T(sim_method.L)
        )
        @views sig[sample, :] .= sens * Mxy[:, i]
    end
    return nothing
end

function acquire_signal!(
    sig,
    sample,
    M,
    ::NoMotion,
    x,
    y,
    z,
    sim_method::BlochSimplewithMultiCoils,
    sens,
)
    sig[sample, :] .= sens * M.xy
end

function acquire_signal!(
    sig,
    sample,
    M,
    ::Union{Motion{T}, MotionList{T}},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    sim_method::BlochSimplewithMultiCoils,
    ::Nothing,
) where {T<:Real}
    sens = birdcage_sensitivities(
        vec(x), vec(y), vec(z);
        ncoils=sim_method.ncoils,
        radius=T(sim_method.radius),
        L=T(sim_method.L)
    )
    sig[sample, :] .= sens * M.xy
end
