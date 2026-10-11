"""Stores preallocated structs for use in Bloch CPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochCPUPrealloc{
    MType<:Mag,
    BzOldType<:AbstractVector,
    BzNewType<:AbstractVector,
    PhiType<:AbstractVector,
    RotType<:Spinor,
    ΔBzType<:AbstractVector,
    S,
    P,
} <: PreallocResult
    M::MType
    Bz_old::BzOldType
    Bz_new::BzNewType
    ϕ::PhiType
    Rot::RotType
    ΔBz::ΔBzType
    sens::S
    coordinates::P
end

view_sens(sensitivities, _) = sensitivities
view_sens(sensitivities::CoilSensitivities, i) = CoilSensitivities(
    @view(sensitivities.values[i, :]),
    sensitivities.interpolators,
)

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.M[i],
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Rot[i],
        p.ΔBz[i],
        view_sens(p.sens, i),
        view_motion_coordinates(p.coordinates, i),
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(
    ::Union{Bloch,BlochDict}, backend::KA.CPU, obj, M,
    max_block_length, _max_adc_samples, _groupsize, sys,
)
    T = eltype(obj.ρ)
    sens = prealloc_sensitivities(sys.receiver, obj)
    ΔBz = similar(obj.ρ)
    @. ΔBz = obj.Δw / T(2π * γ)
    return BlochCPUPrealloc(
        Mag(
            similar(M.xy),
            similar(M.z)
        ),
        zeros(T, size(obj.x)),
        zeros(T, size(obj.x)),
        zeros(T, size(obj.x)),
        Spinor(
            similar(M.xy),
            similar(M.xy)
        ),
        ΔBz,
        sens,
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
    )
end

"""
    run_spin_precession(obj, seq, Xt, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_precession! function in BlochSimpleSimulationMethod.jl 
optimized for the CPU. Uses a loop to step through time instead of allocating a matrix of size 
NSpins x seq.t. The Bz_old, Bz_new, ϕ, and Mxy arrays are pre-allocated in run_sim_time_iter! so 
that they can be re-used from block to block.
"""
function run_spin_precession!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::Bloch,
    groupsize,
    backend::KA.CPU,
    prealloc::PreallocResult
)
    T = eltype(p.ρ)
    #Rename arrays
    Bz_old = prealloc.Bz_old
    Bz_new = prealloc.Bz_new
    ϕ = prealloc.ϕ
    Mxy = prealloc.M.xy
    ΔBz = prealloc.ΔBz
    #Initialize
    fill!(ϕ, zero(T))
    block_time = zero(seq.Δt[1])
    sample = _sample_counter(sig)
    all_adc = size(sig, 1) == length(seq.Δt)
    x, y, z = spin_coordinates!(
        prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[1],
    )
    @. Bz_old = x * seq.Gx[1] + y * seq.Gy[1] + z * seq.Gz[1] + ΔBz
    #Simulation
    @trace track_numbers=false for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i + 1],
        )
        #Effective Field
        @. Bz_new = x * seq.Gx[i + 1] + y * seq.Gy[i + 1] + z * seq.Gz[i + 1] + ΔBz
        #Rotation
        @. ϕ += (Bz_old + Bz_new) * T(-π * γ) * seq.Δt[i]
        block_time += seq.Δt[i]
        #Acquired Signal
        if !isempty(sig)
            adc = all_adc || seq.ADC[i + 1]
            acquire_signal!(sig, sample, adc, M, sys, prealloc, p, seq, i,
                block_time, (x, y, z))
            sample += adc
        end
        #Update simulation state
        Bz_old .= Bz_new
    end
    #Final Spin-State
    @. M.xy = M.xy * exp(-block_time / p.T2) * cis(ϕ)
    @. M.z = M.z * exp(-block_time / p.T1) + p.ρ * (T(1) - exp(-block_time / p.T1))
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

"""
    run_spin_excitation!(obj, seq, sig, M, sim_method, backend, prealloc)

Alternate implementation of the run_spin_excitation! function in BlochSimpleSimulationMethod.jl 
optimized for the CPU. Uses preallocation for all arrays to reduce memory usage.
"""
function run_spin_excitation!(
    p::Phantom,
    seq::DiscreteSequence,
    sig::AbstractArray,
    M::Mag,
    sys,
    sim_method::Union{Bloch,BlochDict},
    groupsize,
    backend::KA.CPU,
    prealloc::BlochCPUPrealloc
)
    T = eltype(p.ρ)
    #Rename arrays
    Bz = prealloc.Bz_old
    B = prealloc.Bz_new
    φ_half = prealloc.ϕ
    α = prealloc.Rot.α
    β = prealloc.Rot.β
    ΔBz = prealloc.ΔBz
    Maux_xy = prealloc.M.xy
    Maux_z = prealloc.M.z
    #Initialize
    sample = _sample_counter(sig)
    all_adc = size(sig, 1) == length(seq.Δt)
    # Rotating frame -> RF frame
    ψ_start = seq.ψ[1]
    @trace track_numbers=false if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    #Simulation
    @trace track_numbers=false for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i],
        )
        #Effective field
        B1 = seq.B1[i]
        @. Bz = (seq.Gx[i] * x + seq.Gy[i] * y + seq.Gz[i] * z) + ΔBz - seq.Δf[i] / T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        @. B = sqrt(abs2(B1) + Bz^2)
        #Spinor Rotation
        @. φ_half = T(-π * γ) * (B * seq.Δt[i]) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        @. α = cos(φ_half)
        @. B = ifelse(B == 0, T(-π * γ) * seq.Δt[i], sin(φ_half) / B)
        @. α -= complex(zero(Bz), Bz * B)
        @. β = complex(imag(B1) * B, -real(B1) * B)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset_at!(M, seq.t, i + 1, p.motion; replace_by=p.ρ)
        #Acquire signal
        if !isempty(sig)
            adc = all_adc || seq.ADC[i + 1] # ADC at the end of the time step
            acquire_signal!(sig, sample, adc, M, sys, prealloc, p, seq, i, sim_method)
            sample += adc
        end
    end
    # RF frame -> Rotating frame
    ψ_end = seq.ψ[end]
    @trace track_numbers=false if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end
