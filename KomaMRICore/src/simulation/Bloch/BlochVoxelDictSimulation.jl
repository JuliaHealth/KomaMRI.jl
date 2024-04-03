import KomaMRICore: SimulationMethod, initialize_spins_state, mul!, output_Ndim,
    run_spin_excitation!, run_spin_precession!, sim_output_dim


Base.@kwdef struct BlochVoxelDict <: SimulationMethod 
    N_spins_per_voxel::Int=1
end


function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochVoxelDict) where {T<:Real}
    N_voxels = length(obj) ÷ sim_method.N_spins_per_voxel
    return (sum(seq.ADC.N), N_voxels)
end


function run_spin_precession!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
    M::Mag{T}, sim_method::BlochVoxelDict) where {T<:Real}
    #Simulation
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    #Effective field
    Bz = xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz' .+ p.Δw / T(2π * γ)
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π * γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π * γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy precession and relaxation, and Mz relaxation
    tp = cumsum(seq.Δt) # t' = t - t0
    dur = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy = [M.xy M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2)] #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    M.z .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    
    #Acquired signal
    # get indexes of the voxels being processed by the present thread
    idx = floor.(Int, p.Dλ1)
    # for each voxel in the voxel-array (list of phantoms) assign the corresponding magnetization value (sum over the spins within a voxel) and normalize
    N_voxels = size(sig, 2)
    #=     println(sum(Mxy[findall(idx .== 1), findall(seq.ADC)];dims=1)/ sim_method.N_spins_per_voxel )
    sleep(3) =#
    for voxel in 1:N_voxels
        sig[:, voxel] += transpose(sum(Mxy[findall(idx .== voxel), findall(seq.ADC)];dims=1)) / sim_method.N_spins_per_voxel 
    end

    return nothing
end


function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochVoxelDict) where {T<:Real}
    #Simulation
    for s ∈ seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        xt = p.x .+ p.ux(p.x, p.y, p.z, s.t)
        yt = p.y .+ p.uy(p.x, p.y, p.z, s.t)
        zt = p.z .+ p.uz(p.x, p.y, p.z, s.t)
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
    return nothing
end