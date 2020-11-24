using MRIsim, JuMP, Ipopt
using MRIsim: get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, dur, get_bvalue, write_diff_fwf, delay

## Parameters
Gmax =  20e-3 #21.7 mT/m
Smax = 50 #52.5 mT/m/ms

f = readlines("F:/qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]
@info δ1, rf180, δ2

τ = (δ1 + rf180 + δ2) * 1e-3 # τ/Nt = Δt => Nt = τ/Δt  
N1 = 600 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 
N2 = floor(Int, N1 * δ2 / δ1)
DIF = Sequence([Grad_fun(x -> 1, δ1*1e-3, N1) delay(rf180*1e-3) Grad_fun(x -> 1, δ2*1e-3, N2); 
                Grad_fun(x -> 1, δ1*1e-3, N1) delay(rf180*1e-3) Grad_fun(x -> 1, δ2*1e-3, N2)])
idx180 = N1 + 1
_, N = size(DIF.GR)
B =  get_Bmatrix(DIF)
SR = get_SRmatrix(DIF)
M =  get_Mmatrix(DIF);

## Options
k = 1
maxwell = true
plots = true
# Optimazation    
Mm = M[1:k+1,:]
model = Model(); set_optimizer(model, Ipopt.Optimizer); #set_silent(model)
@variable(model, -Gmax <= x[1:N] <= Gmax, start=1); #max-grads
@objective(model, Max, x'*B*x); #b-value
@constraint(model, rf180,  x[idx180] .== 0); #rf
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
#TO SCANNER
write_diff_fwf(DIF,idx180,Gmax; filename="F:/qte_vectors_input.txt", name="MXM1")
## Plots
if plots
    MRIsim.plot_grads_moments(DIF,title="ODTI, b=$(round(get_bvalue(DIF)*1e-6, digits=2)) s/mm2")
    ## READ
    G, n1, n2 = MRIsim.read_diff_fwf("F:/qte_vectors_input.txt"); #test read
    using Plots
    plotly()
    plot(G[1:n1,1]*Gmax*1e3, linewidth=2, label="fwf1", legend=:outertopright)
    plot!(G[1:n2,4]*Gmax*1e3, linewidth=2, label="fwf2")
    plot!(diff(G[1:n1,1])*Gmax/6.4e-6, linewidth=2, label="SR1", linestyle=:dashdot)
    plot!(diff(G[1:n2,4])*Gmax/6.4e-6, linewidth=2, label="SR2", linestyle=:dashdot)
    plot!(size=(800,400))
end
