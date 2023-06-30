Base.@kwdef struct BlochMov <: SimulationMethod 
    save_Mz::Bool=false
end

export BlochMov
Base.show(io::IO, s::BlochMov) = begin
	print(io, "BlochMov(save_Mz=$(s.save_Mz))")
end


output_Ndim(sim_method::BlochMov) = 2 #time-points x coils

function get_displacements(p::Phantom{T}, t::AbstractVector{T})where {T<:Real}
    Ns = length(p.x)
    dur = p.dur
    limits = get_pieces_limits(dur, p.K)
    times = time_partitioner(t, dur, limits)

    Δx = hcat(zeros(Ns,1),p.Δx,zeros(Ns,1))
    Δy = hcat(zeros(Ns,1),p.Δy,zeros(Ns,1))
    Δz = hcat(zeros(Ns,1),p.Δz,zeros(Ns,1))

    aux_x = zeros(Ns,1)
    aux_y = zeros(Ns,1)
    aux_z = zeros(Ns,1)

    for i in 1:length(times)
        j = k = i
        while j>(length(limits)-1)
            j -= (length(limits)-1)
        end
        while k>(p.K)
            k -= p.K
        end
        α = (times[i] .- limits[j]) ./ (limits[j+1] - limits[j])
        
        aux_x = hcat(aux_x, Δx[:,k+1]*α' + Δx[:,k]*(1 .- α)') 
        aux_y = hcat(aux_y, Δy[:,k+1]*α' + Δy[:,k]*(1 .- α)') 
        aux_z = hcat(aux_z, Δz[:,k+1]*α' + Δz[:,k]*(1 .- α)') 
    end

    Ux = CuArray(aux_x[:,2:end])
    Uy = CuArray(aux_y[:,2:end])
    Uz = CuArray(aux_z[:,2:end])

    Ux,Uy,Uz
end


function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochMov) where {T<:Real}
    return (sum(seq.ADC.N), 1) #Nt x Ncoils, This should consider the coil info from sys
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom{T}, sim_method::BlochMov) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
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
function run_spin_precession!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
    M::Mag{T}, sim_method::BlochMov) where {T<:Real}
    #Simulation

    xt = p.x 
    yt = p.y 
    zt = p.z 

    #Motion 
    """
    TO-DO: improve function get_displacements to use GPU 
    Ux, Uy, Uz = get_displacements(p,seq.t)

    xt = p.x .+ Ux
    yt = p.y .+ Uy
    zt = p.z .+ Uz
    """

    #Effective field
    Bz = Float32.(xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz' .+ p.Δw / T(2π * γ))
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π * γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π * γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy preccesion and relaxation, and Mz relaxation
    tp = cumsum(seq.Δt) # t' = t - t0
    dur = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy = M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2) #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    #Acquired signal
    sig .= transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1)) #<--- TODO: add coil sensitivities
    return nothing
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
function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochMov) where {T<:Real}
    #Simulation
    for s ∈ seq #This iterates over seq, "s = seq[i,:]"
        xt = p.x 
        yt = p.y 
        zt = p.z 

        #Motion
        """ TO-DO: improve function get_displacements to use GPU
        Ux, Uy, Uz = get_displacements(p,s.t)


        if size(Ux)[2] == 0
            Ux = CuArray(zeros(length(p.x)))
        end
        if size(Uy)[2] == 0
            Uy = CuArray(zeros(length(p.y)))
        end
        if size(Uz)[2] == 0
            Uz = CuArray(zeros(length(p.z)))
        end

        xt = p.x + reshape(Ux,(length(p.x),))
        yt = p.y + reshape(Ux,(length(p.y),))
        zt = p.z + reshape(Ux,(length(p.z),))
        """
    
        #Effective field
        ΔBz = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
        Bz = (s.Gx .* xt .+ s.Gy .* yt .+ s.Gz .* zt) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        mul!( Q(φ, s.B1 ./ B, Bz ./ B), M )
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z  .= M.z  .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    #Acquired signal
    #sig .= -0.1im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end


    