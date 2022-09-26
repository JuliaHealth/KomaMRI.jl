# Code used to generate moment-compensated diffusion gradient waveforms
# Sequence optimization for diffusion prepared motion-compensated MRF 

using KomaMRI, JuMP, Ipopt
using LinearAlgebra: I, Bidiagonal, norm
using Printf

## Aux functions
""""Calculates the normalized moments at the end of the sequence. """
function get_Mmatrix(seq::Sequence; axis=1)
    τ = dur(seq)
    T0 = cumsum([0; seq.DUR])
    M0, M1, M2, M3 = Float64[], Float64[], Float64[], Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) #Durations of pulse
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations
            append!(M0, δ/τ)
            append!(M1, δ.*(T .+ δ/2)/τ^2)
            append!(M2, δ.*(T.^2 .+ T.*δ .+ δ.^2/3)/τ^3)
            append!(M3, δ.*(T.^3 .+ 3/2 * T.^2 .*δ .+ T.*δ.^2 .+ δ.^3/4)/τ^4)
        end
    end
    [M0'; M1'; M2'; M3']
end

"""Slew rate matrix: |SR*g| ≤ Smax."""
function get_SRmatrix(seq::Sequence; axis = 1)
    SR = Bidiagonal{Float64}[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        if N > 1
            Δt = ones(N) * Gi.T / (N-1)
            dv = Δt
            ev = Δt[1:end-1]
            SRi = Bidiagonal(-1 ./ dv, 1 ./ ev, :U)
            # SRi = [SRi[1,:]' ; SRi]; SRi[1,1] = 1/Δt[1] 
            push!(SR, SRi)
        end
    end
    SR
end

"Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m]."
get_Bmatrix(seq::Sequence; axis=1) = begin
    T0 = cumsum([0; seq.DUR[1:end-1]])
    #Calculating timings
    T = Float64[]
    δ = Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay
        if N > 1
            δi = ones(N) * Gi.T / (N-1)
            Ti = [sum(δi[1:j]) for j = 1:N-1]
            Ti = T0[i] .+ delay .+ [0; Ti] #Position of pulse
            append!(δ, δi)
            append!(T, Ti)
        end
    end
    τ = dur(seq) + δ[end]
    Nsamples = length(T)
	ij = [max(i,j) for i=1:Nsamples, j=1:Nsamples]
	α = [(i==j) ? 2/3 : 1/2 for i=1:Nsamples, j=1:Nsamples]
	b = (δ' .* δ) .* (τ .- T[ij] .- α .* δ[ij])
	b_value = (2π*γ)^2*1e-6*b # Trace of B tensor
	b_value
end

## TO SCANNER (Philips)
"""Duration in [s] => samples, with dwell-time of Δt = 6.4 μs."""
δ2N(δ) = floor(Int64, δ * 156250) + 2

"""Exports diffusion preparation waveforms for their use in the scanner."""
function write_diffprep_fwf(G1, G2, G3, bmax, Gmax; filename="./qte_vectors_input.txt", name="Maxwell2")
    open(filename, "w") do io
        t1 = range(0, G1.GR[1].T, length=δ2N(G1.GR[1].T))
		t2 = range(0, G2.GR[1].T, length=δ2N(G2.GR[1].T))
        t3 = range(0, G3.GR[1].T, length=δ2N(G3.GR[1].T))
        maxN = max(length(t1), length(t2), length(t3))
        Gx1, Gy1, Gz1 = KomaMRI.get_grads(G1, Array(t1).+G1.GR[1].delay)
		Gx2, Gy2, Gz2 = KomaMRI.get_grads(G2, Array(t2).+G2.GR[1].delay)
        Gx3, Gy3, Gz3 = KomaMRI.get_grads(G3, Array(t3).+G3.GR[1].delay)
        line = @sprintf "%s %i %i %i % .4f % .4f\n" name length(t1) length(t2) length(t3) bmax Gmax*1e3
        write(io, line)
        for i = 1:maxN
            fx1, fy1, fz1 = i ≤ length(t1) ? (Gx1[i], Gy1[i], Gz1[i])./Gmax   : (0,0,0)
            fx2, fy2, fz2 = i ≤ length(t2) ? (Gx2[i], Gy2[i], Gz2[i])./Gmax   : (0,0,0)
            fx3, fy3, fz3 = i ≤ length(t3) ? (Gx3[i], Gy3[i], Gz3[i])./Gmax   : (0,0,0)
            line = @sprintf "% .4f % .4f % .4f % .4f % .4f % .4f % .4f % .4f % .4f\n" fx1 fy1 fz1 fx2 fy2 fz2 fx3 fy3 fz3
            write(io, line)
        end
    end
end

## Parameters
Gmax =  40e-3 # T/m
Smax = 200 # mT/m/ms
N1 = 400 # You can solve the opt problem in a lower time resolution or use δ2N(dur_grad) 
path_file = "/home/ccp/"
k = 2 #Number of moments to null
maxwell = true #maxwell or concomitant gradient compensation
seq_name = maxwell ? "MX_MC$(k)" : "MC$(k)" #Name of the sequnce

# Timings
Δ1, Δ2 = 15.6180e-3, 45.6180e-3
δ1, δ2, δ3 = 13.8004e-3, 28.1824e-3, 13.8004e-3
rf1 = Δ1 - δ1
rf2 = Δ2 - δ2 - Δ1
# Grads - Pre-defined RF waveforms.
τ = Δ2 + δ3 # τ/Nt = Δt => Nt = τ/Δt  
N2 = floor(Int, N1 * δ2 / δ1)
N3 = floor(Int, N1 * δ3 / δ1)
DIF =  Sequence([Grad(x -> 1e-3,   δ1, N1, 0)])
DIF += Sequence([Grad(x -> 1e-3,   δ2, N2, rf1)])
DIF += Sequence([Grad(x -> 1e-3,   δ3, N3, rf2)])
# Opt matrices
B =  get_Bmatrix(DIF)  #B-value
SR = get_SRmatrix(DIF) #Slew-rate matrices
M =  get_Mmatrix(DIF)  #Moments

## Optimazation
Mm = M[1:k+1,:] #./ 10 .^(2:k+1)
model = Model(); set_optimizer(model, Ipopt.Optimizer); set_silent(model)
@variable(model, -Gmax <= g1[1:N1] <= Gmax, start=Gmax); #max-grads
@variable(model, -Gmax <= g2[1:N2] <= Gmax, start=Gmax); #max-grads
@variable(model, -Gmax <= g3[1:N3] <= Gmax, start=Gmax); #max-grads
@objective(model, Max, [g1;g2;g3]'*B*[g1;g2;g3]); #b-value
@constraint(model, moments, Mm*[g1;g2;g3] .== 0); #moments
@constraint(model, slewrate, -Smax .<= [SR[1]*g1; SR[2]*g2; SR[3]*g3] .<= Smax); #slew rate
@constraint(model, ends, [g1[1]; g2[1]; g3[1]; g1[N1]; g2[N2]; g3[N3]] .== 0)
if maxwell
    @constraint(model, concomitant, sum(g1.^2) - sum(g2.^2) + sum(g3.^2) == 0); #concomitant
end
optimize!(model);
gx1 = value.(g1) #retrieving solution
gx2 = value.(g2) #retrieving solution
gx3 = value.(g3) #retrieving solution

## Solution to Sequence object
DIF[1].GR[1].A =  gx1
DIF[2].GR[1].A = -gx2 #This inversion is necessary, the opt. does not consider the rf180
DIF[3].GR[1].A =  gx3
gx = [gx1; gx2; gx3]
bmax = gx'*B*gx
@info "λ0 = $(abs(round(M[1,:]'*gx/Gmax,digits=1))), λ1 = $(abs(round(M[2,:]'*gx/Gmax,digits=1))), λ2 = $(abs(round(M[3,:]'*gx/Gmax,digits=1)))"
@info "b-value: $(round(bmax, digits=2)) s/mm2"
@info seq_name
if termination_status(model) == MOI.LOCALLY_SOLVED
    @info "Solved!" 
else
    @info "NOT Solved :(" 
end
# TO SCANNER
inv = false
DIFinv = inv ? -DIF : DIF
write_diffprep_fwf(DIFinv[1], DIFinv[2], DIFinv[3], bmax, Gmax; filename="/home/ccp/DiffPrepWaveforms/$seq_name.txt", name=seq_name)
# Plots
p = plot_seq(DIFinv; darkmode=false, slider=false)
savefig(p,"/home/ccp/DiffPrepWaveforms/$seq_name.svg")
p
