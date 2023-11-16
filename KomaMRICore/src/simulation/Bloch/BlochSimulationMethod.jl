struct Bloch <: SimulationMethod end

export Bloch
using LinearAlgebra
mutable struct simParametersPhantom
    M0c::Any #[?]
    T1::Any #[ms]
    T2::Any #[ms]
    relax::Any # bool
    gamma::Any # Hz
    B1c::Any # ?
    dt::Any # seconds
    points::Any
    M0::Any # 3xnpoints mat
    B0::Any # T
    B0_var::Any # PPM
end

include("Magnetization.jl") #Defines Mag <: SpinStateRepresentation
@functor Mag #Gives gpu acceleration capabilities, see GPUFunctions.jl
output_Ndim(sim_method::Bloch) = 2 #time-points x coils

function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::Bloch) where {T<:Real}
    return (sum(seq.ADC.N), 1) #Nt x Ncoils, This should consider the coil info from sys
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom{T}, sim_method::Bloch) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
end

# """
#     run_spin_precession(obj, seq, Xt, sig)

# Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
# = ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
# precession.

# # Arguments
# - `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
# - `seq`: (`::Sequence`) Sequence struct

# # Returns
# - `S`: (`Vector{ComplexF64}`) raw signal over time
# - `M0`: (`::Vector{Mag}`) final state of the Mag vector
# """

# function run_spin_precession!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}}, 
#     M::Mag{T}, sim_method::Bloch) where {T<:Real}
#     #Simulation
#     #Motion
#     xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
#     yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
#     zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
#     #Effective field
#     Bz = xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz' .+ p.Δw / T(2π * γ)
#     #Rotation
#     if is_ADC_on(seq)
#         ϕ = T(-2π * γ) .* cumtrapz(seq.Δt', Bz)
#     else
#         ϕ = T(-2π * γ) .* trapz(seq.Δt', Bz)
#     end
#     #Mxy precession and relaxation, and Mz relaxation
#     tp = cumsum(seq.Δt) # t' = t - t0
#     dur = sum(seq.Δt)   # Total length, used for signal relaxation
#     Mxy = [M.xy M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2)] #This assumes Δw and T2 are constant in time
#     M.xy .= Mxy[:, end]
#     M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
#     #Acquired signal
#     sig .= transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1)) #<--- TODO: add coil sensitivities
#     return nothing
# end


# """
#     M0 = run_spin_excitation(obj, seq, M0)

# It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
# (including B1, gradients and off resonance) and with respect to a rotation axis.

# # Arguments
# - `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
# - `seq`: (`::Sequence`) Sequence struct

# # Returns
# - `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (actually, it's
#     a part of the complete Mag vector and it's a part of the initial state for the next
#     precession simulation step)
# """

# function run_spin_excitation!(p::Phantom{T}, seq::DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
#     M::Mag{T}, sim_method::Bloch) where {T<:Real}
#     #Simulation
#     for s ∈ seq #This iterates over seq, "s = seq[i,:]"
#         #Motion
#         xt = p.x .+ p.ux(p.x, p.y, p.z, s.t)
#         yt = p.y .+ p.uy(p.x, p.y, p.z, s.t)
#         zt = p.z .+ p.uz(p.x, p.y, p.z, s.t)
#         #Effective field
#         ΔBz = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
#         Bz = (s.Gx .* xt .+ s.Gy .* yt .+ s.Gz .* zt) .+ ΔBz
#         B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
#         B[B .== 0] .= eps(T)
#         #Spinor Rotation
#         φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
#         mul!( Q(φ, s.B1 ./ B, Bz ./ B), M )
#         #Relaxation
#         M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
#         M.z  .= M.z  .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
#     end
#     #Acquired signal
#     #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
#     return nothing
# end

## LUKE! ##

function bloch_symmetric_splitting!(u, v, w1, m1, p1, d::simParametersPhantom, sig1, Xt, ADC1)
    
    # members of d
    M0c = d.M0c
    T1 = d.T1 .* 1.0e3 #[ms]
    T2 = d.T2 .* 1.0e3 #[ms]
    relax = d.relax
    gamma = d.gamma * 2 * pi *1.0e-6 #[Mrad/s/T]
    B1c = d.B1c .* 1.0e-3
    dt = Array(d.dt) .* 1.0e3
    Np = size(d.points, 2)
    M0 = d.M0
    u = Array(u) .* 1.0e6 #[uT]
    v = Array(v) .* 1.0e6 #[uT]
    w = Array(w1) .* 1.0e3
    m = Array(m1) .* 1.0e3
    p = Array(p1) .* 1.0e3
    ADC = Array(ADC1)
    sig = Array(sig1)
    @assert length(T1) == Np
    @assert length(T2) == Np
    Nu = CUDA.length(u)                                                                # Number of time points 
    Mrot = CUDA.copy(M0)
    Mt = CUDA.copy(M0)
    i = 0             
    # Since length(sig) = num true ADC
    for n ∈ 1:(Nu-1)

        gadt = gamma * dt[n] / 2   # Time step multiplied by gamma
        B1R = B1c * u[n] * gadt # RF pulse with inhomogeneities
        B1I = B1c * (-v[n]) * gadt

        D12 = CUDA.exp.(-1 ./ T2 .* (relax * dt[n]))
        D3 = CUDA.exp.(-1 ./ T1 .* (relax * dt[n]))
        D = CuArray([1 0 0]') .* D12' .+ CuArray([0 1 0]') .* D12' + .+CuArray([0 0 1]') .* D3'
        b = CuArray([0 0 1]') .* (M0c .- M0c * CUDA.exp.(-1 / T1 * relax * dt[n]))                   #right hand    
        K = gadt .* (d.points[1, :] .* w[n] .+ d.points[2, :] .* m[n] .+ d.points[3, :] .* p[n]) # Inital Gs with no inhomogeneities
        phi = -CUDA.sqrt.((B1R^2 + B1I^2) .+ K .^ 2)
        cs = CUDA.cos.(phi)
        si = CUDA.sin.(phi)
        n1 = B1R ./ CUDA.abs.(phi)
        n2 = B1I ./ CUDA.abs.(phi)
        n3 = K ./ CUDA.abs.(phi)
        n1[CUDA.isnan.(n1)] .= 1
        n2[CUDA.isnan.(n2)] .= 0
        n3[CUDA.isnan.(n3)] .= 0

        #rotation matrix, 3x3 

        Bd1 = n1 .* n1 .* (1 .- cs) .+ cs
        Bd2 = n1 .* n2 .* (1 .- cs) .- n3 .* si
        Bd3 = n1 .* n3 .* (1 .- cs) .+ n2 .* si
        Bd4 = n2 .* n1 .* (1 .- cs) .+ n3 .* si
        Bd5 = n2 .* n2 .* (1 .- cs) .+ cs
        Bd6 = n2 .* n3 .* (1 .- cs) .- n1 .* si
        Bd7 = n3 .* n1 .* (1 .- cs) .- n2 .* si
        Bd8 = n3 .* n2 .* (1 .- cs) .+ n1 .* si
        Bd9 = n3 .* n3 .* (1 .- cs) .+ cs


        #Computes the solution to M(r,u) using strang splitting techniques

        Mrot[1, :] = Bd1 .* Mt[1, :] .+ Bd2 .* Mt[2, :] .+ Bd3 .* Mt[3, :]
        Mrot[2, :] = Bd4 .* Mt[1, :] .+ Bd5 .* Mt[2, :] .+ Bd6 .* Mt[3, :]
        Mrot[3, :] = Bd7 .* Mt[1, :] .+ Bd8 .* Mt[2, :] .+ Bd9 .* Mt[3, :]

        # Relaxation

        Mt = D .* Mrot .+ b


        # Second Rotation
        Mrot[1, :] = Bd1 .* Mt[1, :] .+ Bd2 .* Mt[2, :] .+ Bd3 .* Mt[3, :]
        Mrot[2, :] = Bd4 .* Mt[1, :] .+ Bd5 .* Mt[2, :] .+ Bd6 .* Mt[3, :]
        Mrot[3, :] = Bd7 .* Mt[1, :] .+ Bd8 .* Mt[2, :] .+ Bd9 .* Mt[3, :]

        Mt .= Mrot
        if (ADC[n])
            i += 1
            sig[i] = sum(Mrot[1, :]) + sum(Mrot[2, :])im     #transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1))
        end
        #println(round(n / Nu,digits=2))
    end

    Xt.xy .= Mrot[1, :] .+ Mrot[2, :] .* im
    Xt.z .= Mrot[3, :]
    sig1 .= CuArray(sig)
    return nothing
end

#= no gpu
function bloch_symmetric_splitting!(u, v, w, m, p,  d::simParametersPhantom, sig, ADC)
    # members of d
    M0c = d.M0c
    T1 = d.T1
    T2 = d.T2
    relax = d.relax
    gamma = d.gamma
    B1c = d.B1c
    dt = d.dt                                                     
    Np = size(d.points,2)   
    M0 = d.M0                                               
    @assert length(T1) == Np
    @assert length(T2) == Np
    Nu = length(u)                                                                # Number of time points 
    Bd = CuArray{Float64}(undef, (9, Np)) 
    cs = CuArray{Float64}(undef, Np)
    si = CuArray{Float64}(undef, Np)
    n1 = CuArray{Float64}(undef, Np)
    n2 = CuArray{Float64}(undef, Np)
    n3 = CuArray{Float64}(undef, Np)
    phi = CuArray{Float64}(undef, Np)
    Mrot = zeros(Float64, (3, Np))
    Mrot .= copy(M0)
    Mt = M0
    D = fill(zeros(3,3),Np)
    b = fill(zeros(3),Np)
    for n ∈ 1:Nu
        #=
        if (n == Nu)
            gadt = gamma * ndt / 2   # Time step multiplied by gamma
            for m ∈ 1:Npf (block == length(seq)) 
            ndt = seq_block.Δt[length(seq_block)]
        else
            n
                D[m] = diagm([exp(-1 / T2[m] * relax * ndt), exp(-1 / T2[m] * relax * ndt),        #relaxation
                exp(-1 / T1[m] * relax * ndt)])
                b[m] = [0; 0; M0c] - [0; 0; M0c * exp(-1 / T1[m] * relax * ndt)]                     #right hand
            end    
        else
            gadt = gamma * dt[n] / 2   # Time step multiplied by gamma
            for m ∈ 1:Np
                D[m] = diagm([exp(-1 / T2[m] * relax * dt[n]), exp(-1 / T2[m] * relax * dt[n]),        #relaxation
                exp(-1 / T1[m] * relax * dt[n])])
                b[m] = [0; 0; M0c] - [0; 0; M0c * exp(-1 / T1[m] * relax * dt[n])]                     #right hand
            end  
        end    
        =#
        gadt = gamma * dt[n] / 2   # Time step multiplied by gamma
        for m ∈ 1:Np
            D[m] = diagm([exp(-1 / T2[m] * relax * dt[n]), exp(-1 / T2[m] * relax * dt[n]),        #relaxation
            exp(-1 / T1[m] * relax * dt[n])])
            b[m] = [0; 0; M0c] - [0; 0; M0c * exp(-1 / T1[m] * relax * dt[n])]                     #right hand
        end                                                               
        B1R = B1c * u[n] * gadt                                                     # RF pulse with inhomogeneities
        B1I = B1c * (-v[n]) * gadt
        K = gadt * (d.points[1,:] .* w[n] + d.points[2,:] .* m[n] + d.points[3,:] .* p[n]) # Inital Gs with no inhomogeneities
        for m ∈ 1:Np
            phi[m] = -sqrt(B1R^2 + B1I^2 + K[m]^2)
        end
        for m ∈ 1:Np
            cs[m] = cos(phi[m])
            si[m] = sin(phi[m])
            n1[m] = B1R / abs(phi[m])
            n2[m] = B1I / abs(phi[m])
            n3[m] = K[m] / abs(phi[m])
        end 
        replace!(n1, NaN => 1)
        replace!(n2, NaN => 0)
        replace!(n3, NaN => 0) 

        #rotation matrix, 3x3 

        for m ∈ 1:Np
            Bd[1,m] = n1[m] * n1[m] * (1 - cs[m]) + cs[m]                       
            Bd[2,m] = n1[m] * n2[m] * (1 - cs[m]) - n3[m] * si[m]
            Bd[3,m] = n1[m] * n3[m] * (1 - cs[m]) + n2[m] * si[m]
            Bd[4,m] = n2[m] * n1[m] * (1 - cs[m]) + n3[m] * si[m]
            Bd[5,m] = n2[m] * n2[m] * (1 - cs[m]) + cs[m]
            Bd[6,m] = n2[m] * n3[m] * (1 - cs[m]) - n1[m] * si[m]
            Bd[7,m] = n3[m] * n1[m] * (1 - cs[m]) - n2[m] * si[m]
            Bd[8,m] = n3[m] * n2[m] * (1 - cs[m]) + n1[m] * si[m]
            Bd[9,m] = n3[m] * n3[m] * (1 - cs[m]) + cs[m]
        end

        #Computes the solution to M(r,u) using strang splitting techniques

        for x ∈ 1:Np
            Mrot[1, x] = Bd[1,x] * Mt[1, x] + Bd[2,x] * Mt[2, x] + Bd[3,x] * Mt[3, x]
            Mrot[2, x] = Bd[4,x] * Mt[1, x] + Bd[5,x] * Mt[2, x] + Bd[6,x] * Mt[3, x]
            Mrot[3, x] = Bd[7,x] * Mt[1, x] + Bd[8,x] * Mt[2, x] + Bd[9,x] * Mt[3, x]
        end

        #Idk why but turbo throws error here
        for x ∈ 1:Np
            Mt[1, x] = D[x][1, 1] * Mrot[1, x]
            Mt[2, x] = D[x][2, 2] * Mrot[2, x]
            Mt[3, x] = D[x][3, 3] * Mrot[3, x] + b[x][3]
        end

        # Second Rotation
        for x ∈ 1:Np
            Mrot[1, x] = Bd[1,x] * Mt[1, x] + Bd[2,x] * Mt[2, x] + Bd[3,x] * Mt[3, x]
            Mrot[2, x] = Bd[4,x] * Mt[1, x] + Bd[5,x] * Mt[2, x] + Bd[6,x] * Mt[3, x]
            Mrot[3, x] = Bd[7,x] * Mt[1, x] + Bd[8,x] * Mt[2, x] + Bd[9,x] * Mt[3, x]
        end
        Mt .= Mrot
        if (ADC[n])
            sig[n] = sum(Mrot[1,:]) + sum(Mrot[2,:])im     #transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1))
        end  
        println(round(n / Nu,digits=2))
    end  
    return nothing
end
=#

function bloch_symmetric_splitting!(u, v, w1, m1, p1, d::simParametersPhantom)

    # Input Vars
    u = Array(u) .* 1.0e6 # [uT]
    v = Array(v) .* 1.0e6 # [uT]
    w = Array(w1) .* 1.0e3 # [mT/m]
    m = Array(m1) .* 1.0e3 # [mT/m]
    p = Array(p1) .* 1.0e3 # [mT/m]

    Np = size(d.points, 2)
    gamma = d.gamma * 2.0 * pi * 1.0e-6 # [Mrad/s]
    
    # members of d
    M0c = d.M0c
    T1 = d.T1 .* 1.0e3 # [ms]
    T2 = d.T2 .* 1.0e3 # [ms]
    relax = d.relax
    
    B1c = d.B1c .* 1e-3
    
    Nu = CUDA.length(u)   
    M0 = d.M0
    dt = Array(d.dt) .* 1.0e3 # [ms]

    @assert length(T1) == Np
    @assert length(T2) == Np
    Mrot = CUDA.copy(M0)
    Mt = CUDA.copy(M0)

    for n ∈ 1:(Nu-1)

        gadt = gamma * dt[n] / 2.0   # Time step multiplied by gamma
        B1R = B1c * u[n] * gadt                                                     # RF pulse with inhomogeneities
        B1I = B1c * v[n] * gadt # need to check if this is + or - from matlab conversion

        D12 = CUDA.exp.(-1.0 ./ T2 .* (relax * dt[n]))
        D3 = CUDA.exp.(-1.0 ./ T1 .* (relax * dt[n]))
        
        D = CuArray([1.0 0.0 0.0]') .* D12' .+ CuArray([0.0 1.0 0.0]') .* D12' + .+CuArray([0.0 0.0 1.0]') .* D3'
        
        b = CuArray([0.0 0.0 1.0]') .* (M0c .- M0c * CUDA.exp.(-1.0 / T1 * relax * dt[n]))                   #right hand    
        
        K = gadt .* (d.points[1, :] .* w[n] .+ d.points[2, :] .* m[n] .+ d.points[3, :] .* p[n]) # Inital Gs with no inhomogeneities
        
        # MATLABic
        phi = -CUDA.sqrt.((B1R.^2 + B1I.^2) .+ K .^ 2)
        cs = CUDA.cos.(phi)
        si = CUDA.sin.(phi)
        
        # MATLABic
        n1 = B1R ./ CUDA.abs.(phi)
        n2 = B1I ./ CUDA.abs.(phi)
        n3 = K ./ CUDA.abs.(phi)
        
        # MATLABic
        n1[CUDA.isnan.(n1)] .= 1.0
        n2[CUDA.isnan.(n2)] .= 0.0
        n3[CUDA.isnan.(n3)] .= 0.0

        #rotation matrix, 3x3 

        Bd1 = n1 .* n1 .* (1.0 .- cs) .+ cs
        Bd2 = n1 .* n2 .* (1.0 .- cs) .- n3 .* si
        Bd3 = n1 .* n3 .* (1.0 .- cs) .+ n2 .* si
        Bd4 = n2 .* n1 .* (1.0 .- cs) .+ n3 .* si
        Bd5 = n2 .* n2 .* (1.0 .- cs) .+ cs
        Bd6 = n2 .* n3 .* (1.0 .- cs) .- n1 .* si
        Bd7 = n3 .* n1 .* (1.0 .- cs) .- n2 .* si
        Bd8 = n3 .* n2 .* (1.0 .- cs) .+ n1 .* si
        Bd9 = n3 .* n3 .* (1.0 .- cs) .+ cs


        #Computes the solution to M(r,u) using strang splitting techniques

        Mrot[1, :] = Bd1 .* Mt[1, :] .+ Bd2 .* Mt[2, :] .+ Bd3 .* Mt[3, :]
        Mrot[2, :] = Bd4 .* Mt[1, :] .+ Bd5 .* Mt[2, :] .+ Bd6 .* Mt[3, :]
        Mrot[3, :] = Bd7 .* Mt[1, :] .+ Bd8 .* Mt[2, :] .+ Bd9 .* Mt[3, :]


        Mt = D .* Mrot .+ b


        # Second Rotation
        Mrot[1, :] = Bd1 .* Mt[1, :] .+ Bd2 .* Mt[2, :] .+ Bd3 .* Mt[3, :]
        Mrot[2, :] = Bd4 .* Mt[1, :] .+ Bd5 .* Mt[2, :] .+ Bd6 .* Mt[3, :]
        Mrot[3, :] = Bd7 .* Mt[1, :] .+ Bd8 .* Mt[2, :] .+ Bd9 .* Mt[3, :]

        Mt .= Mrot
    end
    return Mt
end
