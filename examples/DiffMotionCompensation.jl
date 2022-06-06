# Code used to generate gradient waveforms that are moment-compensated
# Sequence optimization for diffusion motion-compensation 

using KomaMRI, JuMP, Ipopt #, Gtk
using KomaMRI: get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              dur, get_bvalue, write_diff_fwf, delay
## Parameters
dwell_time = 6.4e-6
Gmax =  21e-3 #21e-3 45e-3 T/m
Smax = 50 #52.5 mT/m/ms
plots = true 
N1 = 250 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 
#Path were to write the waveforms
#path = "/media/ccp/Samsung_T5/"
path_file = "/home/ccp/Documents/KomaMRI.jl/"
k = 2
sym = false
maxwell = true
moment_nulled, bvalue0 = false, 50
seq_name = maxwell ? "MX_MC$(k)_" : "MC$(k)_"
seq_name = sym ? seq_name*"sym_" : seq_name
seq_name = seq_name * "b" * string(bvalue0)
# Read files
f = readlines(path_file*"qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]*1e-3
if sym 
    rf180 += δ1 - δ2
    δ1 = δ2
end
#Grads# - Pre-defined RF waveforms.
τ = (δ1 + rf180 + δ2) # τ/Nt = Δt => Nt = τ/Δt  
N2 = floor(Int, N1 * δ2 / δ1)
DIF = Sequence([Grad(x -> 1, δ1 - dwell_time, N1) Delay(rf180) Grad(x -> 1, δ2 - dwell_time, N2)])
idx180 = N1 + 1
_, N = size(DIF.GR)
#Opt matrices
B =  get_Bmatrix(DIF)  #B-value
SR = get_SRmatrix(DIF) #Slew-rate
M =  get_Mmatrix(DIF; order=0); #Moments
M0v, M1v, M2v = M[1,:], M[2,:], M[3,:]
t = DIF.GR.T
T = [sum(t[1:i])  for i=1:N] #Acumulated time
# Optimazation
Mm = M[2:k+1,:] ./ 10 .^(2:k+1)
println("###")
model = Model(); set_optimizer(model, Ipopt.Optimizer); set_silent(model)
@variable(model, -Gmax <= x[1:N] <= Gmax, start=1, start=Gmax); #max-grads

if moment_nulled
    @objective(model, Min, -x'*B*x); #b-value
    @constraint(model, moments, Mm*x .== 0); #moments
else
    @variable(model, t1)
    # @variable(model, t2)
    @objective(model, Min, t1);
    @constraint(model, bvalue, x'*B*x == bvalue0); #b-value
    @constraint(model, moments, [t1; Mm*x] in MOI.NormOneCone(k+1))
    # @constraint(model, [t2; 1e-7*SR*x] in MOI.NormOneCone(N+2))
end
@constraint(model, M0, M0v'*x == 0) #diffusion condition
@constraint(model, start_seq,  x[1] .== 0); #seq
@constraint(model, end_seq,  x[end] .== 0); #seq
@constraint(model, start_RF180,  x[idx180] .== 0); #rf
@constraint(model, end_RF180,  x[idx180 + 1] .== 0); #rf
@constraint(model, slewrate, -Smax .<= SR*x .<= Smax); #slew rate
if maxwell
    @constraint(model, concomitant, sum(x[1:idx180-1].^2) == sum(x[idx180+1:end].^2)); #concomitant
end
optimize!(model);
gx = value.(x) #retrieving solution
#Seq object
setproperty!.(DIF.GR, :A, [gx'; 0*gx'; 0*gx'])
bmax = gx'*B*gx
@info "λ0 = $(abs(round(M0v'*gx/Gmax,digits=1))), λ1 = $(abs(round(M1v'*gx/Gmax,digits=1))), λ2 = $(abs(round(M2v'*gx/Gmax,digits=1)))"
@info "b-value: $(round(bmax, digits=2)) s/mm2"
@info round(δ1*1e3,digits=4), round(rf180*1e3,digits=4), round(δ2*1e3,digits=4)
@info seq_name
if termination_status(model) == MOI.LOCALLY_SOLVED
    @info "Solved!" 
else
    @info "NOT Solved :(" 
end
#TO SCANNER
write_diff_fwf(DIF,idx180,Gmax,floor(Int64,bmax); filename=path_file*"QTE_Waveforms/$(seq_name).txt", name=seq_name)
write_diff_fwf(DIF,idx180,Gmax,floor(Int64,bmax); filename=path_file*"QTE_Waveforms/qte_vectors_input.txt", name=seq_name)
# Plots
# if plots
KomaMRI.plot_grads_moments(DIF,title="ODTI, b=$(round(get_bvalue(DIF)*1e-6, digits=2)) s/mm2, λ0 = $(abs(round(1e3*M0v'*gx,digits=1))), λ1 = $(abs(round(1e3*M1v'*gx,digits=1))), λ2 = $(abs(round(1e3*M2v'*gx,digits=1)))")
# end
# end

