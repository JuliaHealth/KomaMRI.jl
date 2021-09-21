# Code used to generate gradient waveforms that are moment-compensated
# Sequence optimization for diffusion motion-compensation 

using LinearAlgebra: norm1
using MRIsim, JuMP, Ipopt #, Gtk
using MRIsim: get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, dur, get_bvalue, write_diff_fwf, delay
using LinearAlgebra
using MathOptInterface: NormOneCone, NormInfinityCone
## Parameters
dwell_time = 6.4e-6
Gmax =  21e-3 #21e-3 45e-3
Smax = 50 #52.5 mT/m/ms
plots = true 
N1 = 400 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 

#Path were to write the waveforms
path = "/home/ccp/Documents/MRIsim.jl/"
#path = "/media/ccp/9549-1842/"

k = 2
sym = false
maxwell = true
# for bvalue0 = [795] #[200, 400, 600, 800]
#for k = [0,1,2], (sym,maxwell)=[(true,false),(false,false),(false,true)]
seq_name = maxwell ? "MXM$k" : "M$k"
seq_name = sym ? seq_name*"_sym" : seq_name
# seq_name = seq_name * string(bvalue0)
# Read files
f = readlines(path*"qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]*1e-3
if sym 
    rf180 += δ1 - δ2
    δ1 = δ2
end
#Grads# - Pre-defined RF waveforms.

τ = (δ1 + rf180 + δ2) # τ/Nt = Δt => Nt = τ/Δt  
N2 = floor(Int, N1 * δ2 / δ1)
DIF = Sequence([Grad_fun(x -> 1, δ1 - dwell_time, N1) delay(rf180) Grad_fun(x -> 1, δ2 - dwell_time, N2); 
                Grad_fun(x -> 1, δ1 - dwell_time, N1) delay(rf180) Grad_fun(x -> 1, δ2 - dwell_time, N2)])
idx180 = N1 + 1
_, N = size(DIF.GR)
#Opt matrices
B =  get_Bmatrix(DIF)  #B-value
SR = get_SRmatrix(DIF) #Slew-rate
M =  get_Mmatrix(DIF); #Moments
M0v = M[1,:]
M1v = M[2,:]
M2v = M[3,:]
# M12v = M[2:3,:]
M0 = [j<=i ? M0v[j] : 0 for i=1:N, j=1:N]
t = DIF.GR.T
T = [sum(t[1:i])  for i=1:N] #Acumulated time
# Optimazation
# λ = 5
Mm = M[1:k+1,:]
println("###")
model = Model(); set_optimizer(model, Ipopt.Optimizer); set_silent(model)
@variable(model, -Gmax <= x[1:N] <= Gmax, start=1, start=Gmax); #max-grads
# @variable(model, t)
@objective(model, Min, -x'*B*x); #+); #b-value
# @objective(model, Min, t); #+); #b-value
# @constraint(model, [t; M1v'*x] in NormOneCone(2))
# @constraint(model, [t; 1e5*1e-6*SR*x] in NormOneCone(N+2))
@constraint(model, start_seq,  x[1] .== 0); #seq
@constraint(model, end_seq,  x[end] .== 0); #seq
@constraint(model, start_RF180,  x[idx180] .== 0); #rf
@constraint(model, end_RF180,  x[idx180 + 1] .== 0); #rf
@constraint(model, slewrate, -Smax .<= SR*x .<= Smax); #slew rate
@constraint(model, moments, Mm*x .== 0); #moments
# @constraint(model, moments1, M1v'*x .== 0); #moments
# @constraint(model, moments2, M2v'*x .== 0); #moments
# @constraint(model, M1extra, -λ <= 1e3*M1v'*x <= λ); #moments
# @constraint(model, bvalue, x'*B*x == bvalue0); #moments
if maxwell
    @constraint(model, concomitant, sum(x[1:idx180-1].^2) == sum(x[idx180+1:end].^2)); #concomitant
end
optimize!(model);
gx = value.(x) #retrieving solution
#Seq object
bmax = gx'*B*gx
setproperty!.(DIF.GR, :A, [gx'; 0*gx'])
@info "λ0 = $(abs(round(1e3*M0v'*gx,digits=1))), λ1 = $(abs(round(1e3*M1v'*gx,digits=1))), λ2 = $(abs(round(1e3*M2v'*gx,digits=1)))"
@info "b-value: $(round(bmax, digits=2)) s/mm2"
@info round(δ1*1e3,digits=4), round(rf180*1e3,digits=4), round(δ2*1e3,digits=4)
@info seq_name
#TO SCANNER
write_diff_fwf(DIF,idx180,Gmax,floor(Int64,bmax); filename=path*"QTE_Waveforms/$(seq_name).txt", name=seq_name)
write_diff_fwf(DIF,idx180,Gmax,floor(Int64,bmax); filename=path*"QTE_Waveforms/qte_vectors_input.txt", name=seq_name)
# Plots
# if plots
MRIsim.plot_grads_moments(DIF,title="ODTI, b=$(round(get_bvalue(DIF)*1e-6, digits=2)) s/mm2, λ = $(abs(round(1e3*M1v'*gx,digits=1)))")
# end
# end

