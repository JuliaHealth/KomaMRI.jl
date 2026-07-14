include("prealloc/Common.jl")
include("prealloc/MagnusConstPrealloc.jl")
include("prealloc/MagnusLinPrealloc.jl")
include("prealloc/MagnusMidPrealloc.jl")
include("prealloc/MagnusQuadPrealloc.jl")
include("prealloc/MagnusGLPrealloc.jl")
include("prealloc/MagnusBGLPrealloc.jl")

# Use Bloch implementation for precession
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnus,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusCPUPrealloc{T}
) where {T<:Real}
    Bz_0, Bz_1 = precession_buffers(prealloc)
    ϕ = prealloc.rotation_norm
    Mxy = prealloc.Maux_xy
    ΔBz = prealloc.ΔBz

    fill!(ϕ, zero(T))
    block_time = zero(T)
    sample = 1
    x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[1])
    @. Bz_0 = x * seq.Gx[1] + y * seq.Gy[1] + z * seq.Gz[1] + ΔBz
    for i in eachindex(seq.Δt)
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        @. Bz_1 = x * seq.Gx[i + 1] + y * seq.Gy[i + 1] + z * seq.Gz[i + 1] + ΔBz
        @. ϕ += (Bz_0 + Bz_1) * T(-π * γ) * seq.Δt[i]
        block_time += seq.Δt[i]
        if seq.ADC[i + 1]
            @. Mxy = exp(-block_time / p.T2) * M.xy * cis(ϕ)
            outflow_spin_reset!(Mxy, seq.t[i + 1], p.motion)
            acquire_signal!(@view(sig[sample, :]), p, sys.receiver, Mxy)
            sample += 1
        end
        Bz_0 .= Bz_1
    end
    @. M.xy = M.xy * exp(-block_time / p.T2) * cis(ϕ)
    @. M.z = M.z * exp(-block_time / p.T1) + p.ρ * (T(1) - exp(-block_time / p.T1))
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ)
    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::BlochMagnusConst1,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusConstCPUPrealloc{T}
) where {T<:Real}
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_0, ωz_0, ωz_1, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    @. ωxy_0 *= B_to_ω
    @. ωz_0  *= B_to_ω
    #Initialize
    sample = 1
    # Rotating frame -> RF frame
    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[1])
    @. ωxy_0 = seq.B1[1] * B_to_ω
    @. ωz_0 = (seq.Gx[1] * x + seq.Gy[1] * y + seq.Gz[1] * z + ΔBz) * B_to_ω + seq.Δf[1] * T(2π)
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        rotation_vector!(θxy, θz, ωxy_0, ωz_0, seq.Δt[i], sim_method)
        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.
        #Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset_at!(M, seq.t, i + 1, p.motion; replace_by=p.ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            acquire_signal!(@view(sig[sample, :]), p, sys.receiver, M.xy)
            sample += 1
        end
        #Update simulation state
        @. ωxy_0 = seq.B1[i + 1] * B_to_ω
        @. ωz_1  = (seq.Gx[i + 1] * x + seq.Gy[i + 1] * y + seq.Gz[i + 1] * z + ΔBz) * B_to_ω + seq.Δf[i + 1] * T(2π)
        ωz_0 .= ωz_1
    end
    @. prealloc.ωxy_0 = ωxy_0 / B_to_ω
    @. prealloc.ωz_0  = ωz_0  / B_to_ω
    # RF frame -> Rotating frame
    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sys,
    sim_method::Union{BlochMagnusLin2,BlochMagnusLinComm2},
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusLinCPUPrealloc{T}
) where {T<:Real}
    B_to_ω = T(-2π * γ)
    ΔBz = prealloc.ΔBz
    (; ωxy_0, ωz_0, ωxy_1, ωz_1, θxy, θz, rotation_norm, α, β, Maux_xy, Maux_z) = prealloc
    @. ωxy_0 *= B_to_ω
    @. ωz_0  *= B_to_ω
    #Initialize
    sample = 1
    # Rotating frame -> RF frame
    ψ_start = seq.ψ[1]
    if !iszero(ψ_start)
        @. M.xy = M.xy * cis(-ψ_start)
    end
    x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[1])
    @. ωxy_0 = seq.B1[1] * B_to_ω
    @. ωz_0 = (seq.Gx[1] * x + seq.Gy[1] * y + seq.Gz[1] * z + ΔBz) * B_to_ω + seq.Δf[1] * T(2π)
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = spin_coordinates(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        #Effective field
        @. ωxy_1 = seq.B1[i + 1] * B_to_ω
        @. ωz_1  = (seq.Gx[i + 1] * x + seq.Gy[i + 1] * y + seq.Gz[i + 1] * z + ΔBz) * B_to_ω + seq.Δf[i + 1] * T(2π)
        rotation_vector!(θxy, θz, ωxy_0, ωz_0, ωxy_1, ωz_1, seq.Δt[i], sim_method)
        #Spinor Rotation
        set_rotation_spinor!(α, β, θxy, θz)
        calc_mag_norm!(rotation_norm, M)
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        restore_mag_norm!(rotation_norm, M) # For reduced float precision only.
        #Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset_at!(M, seq.t, i + 1, p.motion; replace_by=p.ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            acquire_signal!(@view(sig[sample, :]), p, sys.receiver, M.xy)
            sample += 1
        end
        #Update simulation state
        ωz_0  .= ωz_1
        ωxy_0 .= ωxy_1
    end
    @. prealloc.ωxy_0 = ωxy_0 / B_to_ω
    @. prealloc.ωz_0  = ωz_0  / B_to_ω
    # RF frame -> Rotating frame
    ψ_end = seq.ψ[end]
    if !iszero(ψ_end)
        @. M.xy = M.xy * cis(ψ_end)
    end
    return nothing
end

include("MagnusMidCPU.jl")
include("MagnusQuadCPU.jl")
include("MagnusGLCPU.jl")
include("MagnusBGLCPU.jl")
