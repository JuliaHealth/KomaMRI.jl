"""Stores preallocated structs for use in Bloch CPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochCPUPrealloc{
    MType<:Mag,
    BzOldType<:AbstractVector,
    BzNewType<:AbstractVector,
    PhiType<:AbstractVector,
    RotType<:Spinor,
    ΔBzType<:AbstractVector,
} <: PreallocResult
    M::MType
    Bz_old::BzOldType
    Bz_new::BzNewType
    ϕ::PhiType
    Rot::RotType
    ΔBz::ΔBzType
end

Base.view(p::BlochCPUPrealloc, i::UnitRange) = begin
    @views BlochCPUPrealloc(
        p.M[i],
        p.Bz_old[i],
        p.Bz_new[i],
        p.ϕ[i],
        p.Rot[i],
        p.ΔBz[i]
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(sim_method::Bloch, backend::KA.CPU, obj::Phantom, M::Mag, max_block_length::Integer, groupsize)
    T = eltype(obj.ρ)
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
        ΔBz
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
    sim_method::Bloch,
    groupsize,
    backend::KA.CPU,
    prealloc::PreallocResult
)
    Base.@nospecialize p
    return run_bloch_precession!(
        p.motion, p.x, p.y, p.z, p.ρ, p.T1, p.T2, seq, sig, M, prealloc,
    )
end

Base.@noinline function run_bloch_precession!(motion, x0, y0, z0, ρ, T1, T2, seq, sig, M, prealloc)
    T = eltype(ρ)
    #Rename arrays
    Bz_old = prealloc.Bz_old
    Bz_new = prealloc.Bz_new
    ϕ = prealloc.ϕ
    Mxy = prealloc.M.xy
    ΔBz = prealloc.ΔBz
    phase_scale = T(-π * γ)
    #Initialize
    fill!(ϕ, zero(T))
    block_time = zero(T)
    sample = 1
    x, y, z = spin_coordinates(motion, x0, y0, z0, seq.t[1])
    set_precession_field!(Bz_old, x, y, z, ΔBz, seq.Gx[1], seq.Gy[1], seq.Gz[1])
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(motion, x0, y0, z0, seq.t[i + 1])
        #Effective Field
        set_precession_field!(
            Bz_new, x, y, z, ΔBz, seq.Gx[i + 1], seq.Gy[i + 1], seq.Gz[i + 1],
        )
        #Rotation
        accumulate_precession_phase!(ϕ, Bz_old, Bz_new, phase_scale, seq.Δt[i])
        block_time += seq.Δt[i]
        #Acquired Signal
        if seq.ADC[i + 1]
            #Update signal
            @. Mxy = exp(-block_time / T2) * M.xy * cis(ϕ)
            #Reset Spin-State (Magnetization). Only for FlowPath
            outflow_spin_reset!(Mxy, seq.t[i + 1], motion)
            #Acquired signal
            sig[sample] = sum(Mxy) 
            sample += 1
        end
        #Update simulation state
        Bz_old, Bz_new = Bz_new, Bz_old
    end
    #Final Spin-State
    @. M.xy = M.xy * exp(-block_time / T2) * cis(ϕ)
    @. M.z = M.z * exp(-block_time / T1) + ρ * (T(1) - exp(-block_time / T1))
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', motion; replace_by=ρ)
    return nothing
end

Base.@noinline function set_precession_field!(Bz, x, y, z, ΔBz, Gx, Gy, Gz)
    @. Bz = x * Gx + y * Gy + z * Gz + ΔBz
    return nothing
end

Base.@noinline function accumulate_precession_phase!(ϕ, Bz_old, Bz_new, scale, Δt)
    @. ϕ += (Bz_old + Bz_new) * scale * Δt
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
    sim_method::Bloch,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochCPUPrealloc
)
    Base.@nospecialize p
    return run_bloch_excitation!(
        p.motion, p.x, p.y, p.z, p.ρ, p.T1, p.T2, seq, sig, M, prealloc,
    )
end

Base.@noinline function run_bloch_excitation!(motion, x0, y0, z0, ρ, T1, T2, seq, sig, M, prealloc)
    T = eltype(ρ)
    #Rename arrays
    Bz = prealloc.Bz_old
    B = prealloc.Bz_new
    φ_half = prealloc.ϕ
    α = prealloc.Rot.α
    β = prealloc.Rot.β
    ΔBz = prealloc.ΔBz
    Maux_xy = prealloc.M.xy
    Maux_z = prealloc.M.z
    γT = T(γ)
    rotation_scale = T(-π * γ)
    epsT = eps(T)
    oneT = T(1)
    #Initialize
    sample = 1
    # Rotating frame -> RF frame
    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(motion, x0, y0, z0, seq.t[i])
        #Effective field
        B1 = seq.B1[i]
        set_effective_field!(
            Bz,
            B,
            x,
            y,
            z,
            ΔBz,
            seq.Gx[i],
            seq.Gy[i],
            seq.Gz[i],
            seq.Δf[i],
            B1,
            γT,
        )
        #Spinor Rotation
        set_bloch_rotation!(φ_half, α, β, B, Bz, B1, seq.Δt[i], rotation_scale, epsT)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        relax_spins!(M, seq.Δt[i], T1, T2, ρ, oneT)
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset_at!(M, seq.t, i + 1, motion; replace_by=ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            sig[sample] = sum(M.xy) 
            sample += 1
        end
    end
    # RF frame -> Rotating frame
    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end

Base.@noinline function set_effective_field!(Bz, B, x, y, z, ΔBz, Gx, Gy, Gz, Δf, B1, γ)
    @. Bz = Gx * x + Gy * y + Gz * z + ΔBz - Δf / γ
    @. B = sqrt(abs2(B1) + Bz^2)
    return nothing
end

Base.@noinline function set_bloch_rotation!(φ_half, α, β, B, Bz, B1, Δt, scale, ε)
    @. φ_half = scale * (B * Δt)
    @. α = cos(φ_half)
    @. B = sin(φ_half) / (B + (B == 0) * ε)
    @. α -= complex(zero(Bz), Bz * B)
    @. β = complex(imag(B1) * B, -real(B1) * B)
    return nothing
end

Base.@noinline function relax_spins!(M, Δt, T1, T2, ρ, oneT)
    @. M.xy = M.xy * exp(-Δt / T2)
    @. M.z = M.z * exp(-Δt / T1) + ρ * (oneT - exp(-Δt / T1))
    return nothing
end
