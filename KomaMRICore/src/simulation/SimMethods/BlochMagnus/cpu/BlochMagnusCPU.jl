"""Stores preallocated structs for use in Bloch CPU run_spin_precession! and run_spin_excitation! functions."""
struct BlochMagnusCPUPrealloc{T} <: PreallocResult{T}
    M::Mag{T}                               # Mag{T}
    Bxy_old::AbstractVector{Complex{T}}     # Vector{T}(Nspins x 1)
    Bz_old::AbstractVector{T}               # Vector{T}(Nspins x 1)
    Bxy_new::AbstractVector{Complex{T}}     # Vector{T}(Nspins x 1)
    Bz_new::AbstractVector{T}               # Vector{T}(Nspins x 1)
    θxy::AbstractVector{Complex{T}}         # Vector{T}(Nspins x 1)
    θz::AbstractVector{T}                   # Vector{T}(Nspins x 1)
    ϕ::AbstractVector{T}                    # Vector{T}(Nspins x 1)
    Rot::Spinor{T}                          # Spinor{T}
    ΔBz::AbstractVector{T}                  # Vector{T}(Nspins x 1)
end

Base.view(p::BlochMagnusCPUPrealloc, i::UnitRange) = begin
    @views BlochMagnusCPUPrealloc(
        p.M[i],
        p.Bxy_old[i],
        p.Bz_old[i],
        p.Bxy_new[i],
        p.Bz_new[i],
        p.θxy[i],
        p.θz[i],
        p.ϕ[i],
        p.Rot[i],
        p.ΔBz[i]
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(sim_method::BlochMagnus, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real}
    return BlochMagnusCPUPrealloc(
        Mag(
            similar(M.xy),
            similar(M.z)
        ),
        zeros(Complex{T}, size(obj.x)),
        zeros(T, size(obj.x)),
        zeros(Complex{T}, size(obj.x)),
        zeros(T, size(obj.x)),
        zeros(Complex{T}, size(obj.x)),
        zeros(T, size(obj.x)),
        zeros(T, size(obj.x)),
        Spinor(
            similar(M.xy),
            similar(M.xy)
        ),
        obj.Δw ./ T(2π .* γ)
    )
end

# Use Bloch implementation for precession
function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochMagnus,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusCPUPrealloc{T}
) where {T<:Real}
    run_spin_precession!(p, seq, sig, M, Bloch(), groupsize, backend, prealloc)
end

# This part changes a bit more
function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochMagnus,
    groupsize,
    backend::KA.CPU,
    prealloc::BlochMagnusCPUPrealloc{T}
) where {T<:Real}
    #Rename arrays, don't mind this part
    B_to_ω = T(-2π * γ)
    ωxy_old = prealloc.Bxy_old; @. ωxy_old *= B_to_ω
    ωz_old  = prealloc.Bz_old;  @. ωz_old  *= B_to_ω
    ωxy_new = prealloc.Bxy_new
    ωz_new  = prealloc.Bz_new
    θxy = prealloc.θxy
    θz  = prealloc.θz
    θ = prealloc.ϕ
    α = prealloc.Rot.α
    β = prealloc.Rot.β
    ΔBz = prealloc.ΔBz
    Maux_xy = prealloc.M.xy
    Maux_z  = prealloc.M.z
    #Initialize
    sample = 1
    #Simulation
    for i in eachindex(seq.Δt)
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t[i + 1])
        #Effective field
        @. ωxy_new = seq.B1[i + 1] * B_to_ω
        @. ωz_new  = (seq.Gx[i + 1] * x + seq.Gy[i + 1] * y + seq.Gz[i + 1] * z + ΔBz) * B_to_ω + seq.Δf[i + 1] * T(2π)
        effective_rotation_vector!(θxy, θz, ωxy_old, ωz_old, ωxy_new, ωz_new, seq.Δt[i], sim_method)
        #Spinor Rotation
        @. θ = sqrt(abs(θxy) ^ 2 + θz ^ 2)
        @. α = cos(θ / T(2)) - Complex{T}(im) * (θz / θ) * sin(θ / T(2))
        @. β = -Complex{T}(im) * (θxy / θ) * sin(θ / T(2))
        mul!(Spinor(α, β), M, Maux_xy, Maux_z)
        #Relaxation
        @. M.xy = M.xy * exp(-seq.Δt[i] / p.T2)
        @. M.z = M.z * exp(-seq.Δt[i] / p.T1) + p.ρ * (T(1) - exp(-seq.Δt[i] / p.T1))
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset!(M,  seq.t[i + 1, :], p.motion; replace_by=p.ρ)
        #Acquire signal
        if seq.ADC[i + 1] # ADC at the end of the time step
            sig[sample] = sum(M.xy) 
            sample += 1
        end
        #Update simulation state
        ωz_old  .= ωz_new
        ωxy_old .= ωxy_new
    end
    @. prealloc.Bxy_old = ωxy_old / B_to_ω
    @. prealloc.Bz_old  = ωz_old  / B_to_ω
    return nothing
end

function effective_rotation_vector!(θxy, θz, ωxy_old, ωz_old, ωxy_new, ωz_new, Δt, sim_method::BlochMagnus1)
    # Ω1_constant
    @. θxy = ωxy_old * Δt
    @. θz  = ωz_old * Δt
    @. θxy[θxy == 0] = eps(eltype(θz)) # It could be that Bxy[told] or Bxy[tnew] is zero
    return nothing
end

function effective_rotation_vector!(θxy, θz, ωxy_old, ωz_old, ωxy_new, ωz_new, Δt, sim_method::BlochMagnus2)
    # Ω1_linear
    @. θxy = (ωxy_old + ωxy_new) * (Δt / 2)
    @. θz  = (ωz_old + ωz_new) * (Δt / 2)
    return nothing
end

function effective_rotation_vector!(θxy, θz, ωxy_old, ωz_old, ωxy_new, ωz_new, Δt, sim_method::BlochMagnus4)
    # Ω1_linear
    effective_rotation_vector!(θxy, θz, ωxy_old, ωz_old, ωxy_new, ωz_new, Δt, BlochMagnus2())
    # Ω2_linear (this is the complex representation of Δt²/12 (ωₙ₊₁ × ωₙ))
    @. θxy .-= im * (ωxy_new * ωz_old - ωxy_old * ωz_new) * (Δt ^ 2 / 12)
    @. θz  .-= imag(conj(ωxy_old) * ωz_new) * (Δt ^ 2 / 12)
    return nothing
end