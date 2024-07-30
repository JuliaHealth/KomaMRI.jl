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
    pre::BlochGPUPrealloc
) where {T<:Real}
    #Simulation
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')
    
    #Initialize arrays
    seq_block = pre.seq_properties[1]
    Bz = @view pre.Bz[:,1:length(seq.t)]
    Tz = @view pre.Tz[:,1:length(seq.t)-1]
    
    #Effective field
    Bz .= x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz' .+ pre.scaled_Δw

    if seq_block.nADC > 0
        #Rotation
        ϕ = @view pre.ϕ[:,1:length(seq.t)-1]
        Mxy = @view pre.Mxy[:,1:seq_block.nADC]
        cumtrapz!(ϕ, Tz, seq.Δt', Bz)
        ϕ_ADC = @view ϕ[:,seq_block.ϕ_indices]

        if seq_block.first_ADC
            Mxy[:,1] .= M.xy
            Mxy[:,2:end] .= M.xy .* exp.(-seq_block.tp_ADC' ./ p.T2) .* _cis.(ϕ_ADC)
        else
            Mxy .= M.xy .* exp.(-seq_block.tp_ADC' ./ p.T2) .* _cis.(ϕ_ADC)
        end

        #Mxy precession and relaxation, and Mz relaxation
        M.z  .= M.z .* exp.(-seq_block.dur ./ p.T1) .+ p.ρ .* (T(1) .- exp.(-seq_block.dur ./ p.T1))
        M.xy .= M.xy .* exp.(-seq_block.dur ./ p.T2) .* _cis.(ϕ[:,end])
        
        #Acquired signal
        sig .= transpose(sum(Mxy; dims=1))
    else
        #Rotation
        ϕ = @view pre.ϕ[:,1]
        trapz!(ϕ, Tz, seq.Δt', Bz)

        #Mxy precession and relaxation, and Mz relaxation
        M.xy .= M.xy .* exp.(-seq_block.dur ./ p.T2) .* _cis.(ϕ)
        M.z  .= M.z .* exp.(-seq_block.dur ./ p.T1) .+ p.ρ .* (T(1) .- exp.(-seq_block.dur ./ p.T1))
    end

    return nothing
end 