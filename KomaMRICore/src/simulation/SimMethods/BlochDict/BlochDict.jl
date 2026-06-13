Base.@kwdef struct BlochDict <: SimulationMethod
    save_Mz::Bool = false
end

export BlochDict
Base.show(io::IO, s::BlochDict) = begin
    print(io, "BlochDict(save_Mz=$(s.save_Mz))")
end

function sim_output_dim(
    obj::SimulationPhantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochDict
) where {T<:Real}
    out_state_dim = sim_method.save_Mz ? 2 : 1
    return (sum(seq.ADC.N), length(obj), out_state_dim)
end

# To fix BlochDict for CPU parallel execution (#204)
function split_sig_per_thread(sig, i, p, sim_method::BlochDict)
    return @view sig[:, p, :, i]
end

"""
    run_spin_precession(obj, seq, Xt, sig)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
precession.

# Arguments
- `obj`: (`::SimulationPhantom`) simulation phantom view
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_spin_precession!(
    p::SimulationPhantom{T},
    seq::DiscreteSequence{T},
    sig::AbstractArray{Complex{T}},
    M::Mag{T},
    sim_method::BlochDict,
    groupsize,
    backend::KA.Backend,
    prealloc::PreallocResult
) where {T<:Real}
    #Motion
    motion = get_motion(p)
    ρ = get_ρ(p)
    T1 = get_T1(p)
    T2 = get_T2(p)
    x, y, z = get_spin_coords(p, seq.t')
    #Effective field
    Bz = x .* seq.Gx' .+ y .* seq.Gy' .+ z .* seq.Gz' .+ get_Δw(p) ./ T(2π .* γ)
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π .* γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π .* γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy precession and relaxation, and Mz relaxation
    tp = cumsum(seq.Δt) # t' = t - t0
    dur = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy = M.xy .* exp.(-tp' ./ T2) .* cis.(ϕ) #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(Mxy, seq.t[2:end]', motion)
    #Acquire signal
    sig[:, :, 1] .= @views transpose(Mxy[:, findall(seq.ADC[2:end])])

    if sim_method.save_Mz
        Mz = M.z .* exp.(-tp' ./ T1) .+ ρ .* (1 .- exp.(-tp' ./ T1)) #Calculate intermediate points
        #Reset Spin-State (Magnetization). Only for FlowPath
        outflow_spin_reset!(Mz, seq.t[2:end]', motion; replace_by=ρ)
        sig[:, :, 2] .= @views transpose(Mz[:, findall(seq.ADC[2:end])]) #Save state to signal
        M.z .= Mz[:, end]
    else
        M.z .= M.z .* exp.(-dur ./ T1) .+ ρ .* (1 .- exp.(-dur ./ T1)) #Jump to the last point
    end
    #Reset Spin-State (Magnetization). Only for FlowPath
    outflow_spin_reset!(M, seq.t[2:end]', motion; replace_by=ρ)
    return nothing
end

# So we can use the same excitation function than BlochSimple
function acquire_signal!(sig, sample, M, sim_method::BlochDict)
    sig[sample, :, 1] .= M.xy
    if sim_method.save_Mz
        sig[sample, :, 2] .= M.z
    end
end
