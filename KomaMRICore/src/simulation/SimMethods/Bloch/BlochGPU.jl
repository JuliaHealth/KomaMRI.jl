include("KernelFunctions.jl")

"""These properties are redundant with seq.Δt and seq.ADC, but it is much faster
to calculate them on the CPU before the simulation is run."""
struct SeqBlockProperties{T<:Real}
    length::Int64
    nADC::Int64
    first_ADC::Bool
    ADC_indices::AbstractVector{Int64}
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
                length(seq_block.t),
                count(seq_block.ADC),
                false,
                Int64[],
                T[],
                zero(T)
            ))
        else
            ADC_indices = findall(seq_block.ADC) .- 1
            if seq_block.ADC[1]
                deleteat!(ADC_indices, 1)
            end
            tp = cumsum(seq_block.Δt)
            push!(seq_properties, SeqBlockProperties(
                length(seq_block.t),
                count(seq_block.ADC),
                seq_block.ADC[1],
                ADC_indices,
                tp[ADC_indices],
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
    B::AbstractMatrix{T}
    φ::AbstractMatrix{T}
    ΔT1::AbstractMatrix{T}
    ΔT2::AbstractMatrix{T}
    Δϕ::AbstractMatrix{T}
    ϕ::AbstractArray{T}
    Mxy::AbstractMatrix{Complex{T}}
    ΔBz::AbstractVector{T}
end

Base.view(p::BlochGPUPrealloc{T}, i::UnitRange) where {T<:Real} = p
function prealloc_block(p::BlochGPUPrealloc{T}, i::Integer) where {T<:Real}
    seq_block = p.seq_properties[i]

    return BlochGPUPrealloc(
        [seq_block],
        view(p.Bz,:,1:seq_block.length),
        view(p.B,:,1:seq_block.length),
        view(p.φ,:,1:seq_block.length-1),
        view(p.ΔT1,:,1:seq_block.length-1),
        view(p.ΔT2,:,1:seq_block.length-1),
        view(p.Δϕ,:,1:seq_block.length-1),
        seq_block.nADC > 0 ? view(p.ϕ,:,1:seq_block.length-1) : view(p.ϕ,:,1),
        view(p.Mxy,:,1:seq_block.nADC),
        p.ΔBz
    )
end

"""Preallocates arrays for use in run_spin_precession! and run_spin_excitation!."""
function prealloc(sim_method::Bloch, backend::KA.GPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, precalc) where {T<:Real}
    if !(precalc isa BlochGPUPrecalc) @error """Sim method Bloch() does not support calling run_sim_time_iter directly. Use method BlochSimple() instead.""" end

    return BlochGPUPrealloc(
        precalc.seq_properties,
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, T, (size(obj.x, 1), max_block_length-1)),
        KA.zeros(backend, Complex{T}, (size(obj.x, 1), max_block_length)),
        obj.Δw ./ T(2π .* γ)
    )
end

@inline function calculate_precession!(Δϕ::AbstractArray{T}, Δt::AbstractArray{T}, Bz::AbstractArray{T}) where {T<:Real}
    Δϕ .= (Bz[:,2:end] .+ Bz[:,1:end-1]) .* Δt .* T(-π .* γ)
end
@inline function apply_precession!(ϕ::AbstractVector{T}, Δϕ::AbstractArray{T}) where {T<:Real}
    ϕ .= sum(Δϕ, dims=2)
end
function apply_precession!(ϕ::AbstractArray{T}, Δϕ::AbstractArray{T}) where {T<:Real}
    cumsum!(ϕ, Δϕ, dims=2)
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
    #Sequence block info
    seq_block = pre.seq_properties[1]
    
    #Effective field
    pre.Bz .= x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz' .+ pre.ΔBz

    #Rotation
    calculate_precession!(pre.Δϕ, seq.Δt', pre.Bz)
    # Reduces Δϕ Nspins x Nt to ϕ Nspins x Nt, if Nadc = 0, to Nspins x 1 
    apply_precession!(pre.ϕ, pre.Δϕ)

    #Acquired signal
    if seq_block.nADC > 0
        ϕ_ADC = @view pre.ϕ[:,seq_block.ADC_indices]
        if seq_block.first_ADC
            pre.Mxy[:,1] .= M.xy
            pre.Mxy[:,2:end] .= M.xy .* exp.(-seq_block.tp_ADC' ./ p.T2) .* _cis.(ϕ_ADC)
            #Reset Spin-State (Magnetization). Only for FlowPath
            outflow_spin_reset!(pre.Mxy, seq_block.tp_ADC', p.motion; seq_t=seq.t, add_t0=true)
        else
            pre.Mxy .= M.xy .* exp.(-seq_block.tp_ADC' ./ p.T2) .* _cis.(ϕ_ADC)
            #Reset Spin-State (Magnetization). Only for FlowPath
            outflow_spin_reset!(pre.Mxy, seq_block.tp_ADC', p.motion; seq_t=seq.t)
        end

        sig .= transpose(sum(pre.Mxy; dims=1))
    end
    
    #Mxy precession and relaxation, and Mz relaxation
    M.z  .= M.z .* exp.(-seq_block.dur ./ p.T1) .+ p.ρ .* (T(1) .- exp.(-seq_block.dur ./ p.T1))
    M.xy .= M.xy .* exp.(-seq_block.dur ./ p.T2) .* _cis.(pre.ϕ[:,end])

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M, seq.t', p.motion; replace_by=p.ρ)

    return nothing
end

function run_spin_excitation!(
    p::Phantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::Bloch,
    backend::KA.GPU,
    pre::BlochGPUPrealloc
) where {T<:Real}
    #Motion
    x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, seq.t')

    #Effective Field
    pre.Bz .= (x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz') .+ pre.ΔBz .- seq.Δf' ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(x,y,z)
    pre.B .= sqrt.(abs.(seq.B1') .^ 2 .+ abs.(pre.Bz) .^ 2)
    
    #Spinor Rotation
    pre.φ .= T(-π .* γ) .* (pre.B[:,1:end-1] .* seq.Δt') # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
    
    #Relaxation
    pre.ΔT1 .= exp.(-seq.Δt' ./ p.T1)
    pre.ΔT2 .= exp.(-seq.Δt' ./ p.T2)

    #Excitation
    apply_excitation!(backend, 512)(M.xy, M.z, pre.φ, seq.B1, pre.Bz, pre.B, pre.ΔT1, pre.ΔT2, p.ρ, ndrange=size(M.xy,1))
    KA.synchronize(backend)

    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M,  seq.t', p.motion; replace_by=p.ρ) # TODO: reset state inside kernel

    return nothing
end