# Sequence optimization for diffusion motion compensation 

using MRIsim, JuMP, Ipopt, SCS
using MRIsim: get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, dur, get_bvalue, write_diff_fwf, delay

## Parameters
Gmax =  21e-3 #21.7 mT/m
Smax = 50 #52.5 mT/m/ms
k = 2
maxwell = true
seq_name = maxwell ? "MXM$k" : "M$k"
plots = true
N1 = 400 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 

## Read files
#path = "/home/ccp/Documents/MRIsim.jl/"
path = "/media/ccp/9549-1842/"

f = readlines(path*"qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]
@info δ1, rf180, δ2
#Grads
τ = (δ1 + rf180 + δ2) * 1e-3 # τ/Nt = Δt => Nt = τ/Δt  
dwell_time = 6.4e-6
N2 = floor(Int, N1 * δ2 / δ1)
DIF = Sequence([Grad_fun(x -> 1, δ1*1e-3 - dwell_time, N1) delay(rf180*1e-3) Grad_fun(x -> 1, δ2*1e-3 - dwell_time, N2); 
                Grad_fun(x -> 1, δ1*1e-3 - dwell_time, N1) delay(rf180*1e-3) Grad_fun(x -> 1, δ2*1e-3 - dwell_time, N2)])
idx180 = N1 + 1
_, N = size(DIF.GR)
#Opt matrices
B =  get_Bmatrix(DIF)
SR = get_SRmatrix(DIF)
M =  get_Mmatrix(DIF);
# Optimazation    
Mm = M[1:k+1,:]
model = Model(); set_optimizer(model, Ipopt.Optimizer); #set_silent(model)
@variable(model, -Gmax <= x[1:N] <= Gmax, start=1, start=Gmax); #max-grads
@objective(model, Max, x'*B*x); #b-value
@constraint(model, RF180,  x[idx180] .== 0); #rf
@constraint(model, start_seq,  x[1] .== 0); #rf
@constraint(model, end_seq,  x[idx180 + 1] .== 0); #rf
@constraint(model, slewrate, -Smax .<= SR*x .<= Smax); #slewrate
@constraint(model, moments, Mm*x .== 0); #moments
if maxwell
    @constraint(model, concomitant, sum(x[1:idx180-1].^2) == sum(x[idx180+1:end].^2)); #concomitant
end
optimize!(model);
gx = value.(x) #retrieving solution
#seq object
setproperty!.(DIF.GR, :A, [gx'; 0*gx']);
@info "$k-order motion-compensated"
@info "b-value: $(round(gx'*B*gx, digits=2)) s/mm2"
@info δ1, rf180, δ2
@info seq_name
#TO SCANNER
write_diff_fwf(DIF,idx180,Gmax; filename=path*"qte_vectors_input.txt", name=seq_name)

## Plots
if plots
    MRIsim.plot_grads_moments(DIF,title="ODTI, b=$(round(get_bvalue(DIF)*1e-6, digits=2)) s/mm2")
end
