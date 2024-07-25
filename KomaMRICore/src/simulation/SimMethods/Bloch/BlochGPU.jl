"""These properties are redundant with seq.Δt and seq.ADC, but it is much faster
to calculate them on the CPU before the simulation is run."""
struct SeqBlockProperties{T<:Real}
    nADC::Int64
    first_ADC::Bool
    ϕ_indices::AbstractVector{Int64}
    tp_ADC::AbstractVector{T}
    dur::T
end

@functor SeqBlockProperties

"""Stores information added to the prealloc struct for use in run_spin_precession!
and run_spin_excitation!. This information is calculated on the CPU before the
simulation is run."""
struct BlochGPUPrecalc{T} <: PrecalcResult{T}
    seq_properties::AbstractVector{SeqBlockProperties{T}}
end

@functor BlochGPUPrecalc

"""Precalculates sequence properties for use in run_spin_precession!"""
function precalculate(
    sim_method::Bloch, 
    backend::KA.GPU,
    seq::DiscreteSequence{T}, 
    parts::Vector{UnitRange{S}}, 
    excitation_bool::Vector{Bool}
) where {T<:Real,S<:Integer}
    seq_properties = SeqBlockProperties{T}[]
    for (block, p) in enumerate(parts)
        seq_block = @view seq[p]
        if excitation_bool[block]
            push!(seq_properties, SeqBlockProperties(
                0,
                false,
                Int64[],
                T[],
                zero(T)
            ))
        else
            ϕ_indices = findall(seq_block.ADC) .- 1
            if (length(ϕ_indices) > 0 && first(ϕ_indices) == 0)
                deleteat!(ϕ_indices, 1)
            end
            tp = cumsum(seq_block.Δt)
            push!(seq_properties, SeqBlockProperties(
                count(seq_block.ADC),
                seq_block.ADC[1],
                ϕ_indices,
                tp[ϕ_indices],
                last(tp)
            ))
        end
    end

    return BlochGPUPrecalc(seq_properties)
end

"""Stores preallocated structs for use in Bloch GPU run_spin_precession function."""
struct BlochGPUPrealloc{T} <: PreallocResult{T}
    seq_properties::AbstractVector{SeqBlockProperties{T}}
    Bz::AbstractMatrix{T}
    Tz::AbstractMatrix{T}
    ϕ::AbstractMatrix{T}
    Mxy::AbstractMatrix{Complex{T}}
    scaled_Δw::AbstractVector{T}
end

Base.view(p::BlochGPUPrealloc{T}, i::UnitRange) where {T<:Real} = p
function prealloc_block(p::BlochGPUPrealloc{T}, i::Integer) where {T<:Real}
    return BlochGPUPrealloc(
        p.seq_properties[i,:],
        p.Bz,
        p.Tz,
        p.ϕ,
        p.Mxy,
        p.scaled_Δw
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(sim_method::Bloch, backend::KA.GPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, precalc) where {T<:Real}
    if !(precalc isa BlochGPUPrecalc) @error """Sim method Bloch() does not support calling run_sim_time_iter directly. Use method BlochSimple() instead.""" end

    return BlochGPUPrealloc(
        precalc.seq_properties,
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, Complex{T}, (size(obj.x, 1), max_block_length)),
        obj.Δw ./ T(2π .* γ)
    )
end

function run_spin_precession!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.GPU,
    prealloc::BlochGPUPrealloc
) where {T<:Real}
    #Simulation
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    
    #Initialize arrays
    seq_block = prealloc.seq_properties[1]
    dur  = seq_block.dur   # Total length, used for signal relaxation
    ϕ_indices = seq_block.ϕ_indices # Indices of ϕ corresponding to ADC point
    tp = seq_block.tp_ADC
    Bz = @view prealloc.Bz[:,1:length(seq.t)]
    Tz = @view prealloc.Tz[:,1:length(seq.t)-1]
    scaled_Δw = prealloc.scaled_Δw
    #Effective field
    Bz .= x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz' .+ scaled_Δw

    if seq_block.nADC > 0
        #Use preallocated arrays
        ϕ = @view prealloc.ϕ[:,1:length(seq.t)-1]
        Mxy = @view prealloc.Mxy[:,1:seq_block.nADC]
        #Rotation
        cumtrapz!(seq.Δt', Bz, Tz, ϕ)
        ϕ .= ϕ .* T(-2π .* γ)
        if seq_block.first_ADC
            Mxy[:,1] .= M.xy
            Mxy[:,2:end] .= M.xy .* exp.(-tp' ./ p.T2) .* _cis.(ϕ[:,ϕ_indices])
        else
            Mxy .= M.xy .* exp.(-tp' ./ p.T2) .* _cis.(ϕ[:,ϕ_indices])
        end
        #Mxy precession and relaxation, and Mz relaxation
        M.xy .= M.xy .* exp.(-dur ./ p.T2) .* _cis.(ϕ[:,end])
        M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (T(1) .- exp.(-dur ./ p.T1))
        sig .= transpose(sum(Mxy, dims=1))
    else
        #Rotation
        ϕ = T(-2π .* γ) .* trapz!(seq.Δt', Bz, Tz)
        #Mxy precession and relaxation, and Mz relaxation
        M.xy .= M.xy .* exp.(-dur ./ p.T2) .* _cis.(ϕ)
        M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (T(1) .- exp.(-dur ./ p.T1))
    end

    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.GPU,
    prealloc::BlochGPUPrecalc
) where {T<:Real}
    #Simulation
    for s in seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, s.t)
        #Effective field
        ΔBz = p.Δw ./ T(2π .* γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
        Bz = (s.Gx .* x .+ s.Gy .* y .+ s.Gz .* z) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π .* γ) .* (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        mul!(Q(φ, s.B1 ./ B, Bz ./ B), M)
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z .= M.z .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    #Acquired signal
    #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end