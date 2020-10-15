using MRIsim, JuMP, Ipopt
using MRIsim: plot_grads, get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, RF_fun, dur, plot_grads_moments, get_bvalue, get_grads,
              write_diff_fwf, δ2N, read_diff_fwf
## Parameters
Gmax =  20e-3 #21.7 mT/m
Smax = 52 #52.5 mT/m/ms
τ = (18.9+78.7)*1e-3 # τ/Nt = Δt => Nt = τ/Δt  
Nt = 1100 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 
DIF = Sequence([Grad_fun(x -> 1, τ, Nt); 
                Grad_fun(x -> 1, τ, Nt)])
start180T = τ - 18.9e-3 - 8e-3
dur180T = 8e-3
Δt = τ/Nt
start180, end180 = floor(Int64, start180T/Δt), ceil(Int64, (start180T + dur180T)/Δt)
_, N = size(DIF.GR)
B =  get_Bmatrix(DIF)
SR = get_SRmatrix(DIF)
M =  get_Mmatrix(DIF);
## Optimization
k = 2    
Mm = M[1:k+1,:]
model = Model(); set_optimizer(model, Ipopt.Optimizer); #set_silent(model)
@variable(model, -Gmax <= x[1:N] <= Gmax, start=1); #max-grads
@objective(model, Max, x'*B*x); #b-value
@constraint(model, rf180,  x[start180:end180] .== 0); #rf
@constraint(model, start_seq,  x[1] .== 0); #rf
# @constraint(model, end_seq,  x[end] .== 0); #rf
@constraint(model, slewrate, -Smax .<= SR*x .<= Smax); #slewrate
@constraint(model, moments, Mm*x .== 0); #moments
@constraint(model, concomitant, sum(x[1:start180-1].^2) == sum(x[end180+1:end].^2)); #concomitant
optimize!(model);
gx = value.(x) #retrieving solution
#seq object
setproperty!.(DIF.GR, :A, [gx'; 0*gx']);
@info "$k-order motion-compensated"
@info "b-value: $(round(gx'*B*gx, digits=2)) s/mm2"
#Plots
# plot_grads_moments(DIF,title="Maxwell,M0,M1,M2-compensated, b=$(round(get_bvalue(DIF)*1e-6, digits=2)) s/mm2")
plot_grads(DIF)
## TO SCANNER
write_diff_fwf(DIF,start180,end180,Gmax)
G, n1, n2 = read_diff_fwf(); #test read
## READ
using Plots
plotly()
dt = 6.4e-6; #6.4us dwell-time
plot(G[1:n1,1], linewidth=2, label="fwf1", legend=:outertopright)
plot!(G[1:n2,4], linewidth=2, label="fwf2")
plot!(size=(800,400))
ylims!((-1.1,1.1))