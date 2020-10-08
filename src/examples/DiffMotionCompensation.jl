using MRIsim, JuMP, Ipopt
using MRIsim: plot_grads, get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, RF_fun, dur, plot_grads_moments, get_bvalue, get_grads
## Parameters
Gmax =  30e-3 #21.7e-3; 
Smax = 52.5
τ = 70e-3 # τ/Nt = Δt => Nt = τ/Δt  
Nt = 1100 #floor(Int64, τ * 1e3 * 15625 / 100) + 2 # Δt = 6.4e-6 #dwell-time 
DIF = Sequence([Grad_fun(x -> 1, τ, Nt); 
                Grad_fun(x -> 1, τ, Nt)])
start180T = 33e-3
dur180T = 7e-3
Δt = τ/Nt
start180, end180 = floor(Int64, start180T/Δt), floor(Int64, (start180T + dur180T)/Δt)
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
p3 = plot_grads_moments(DIF,title="Maxwell,M0,M1,M2-compensated, b=$(round(get_bvalue(DIF[3])*1e-6, digits=2)) s/mm2")
## TO SCANNER
using Printf
open("./Maxwell2.txt", "w") do io
    dt = 6.4e-6; #6.4us dwell-time
    t = range(0,start180T,step=dt)
    Gx1, _ = get_grads(DIF[1:start180-1],t)
    Gx2, _ = get_grads(DIF[end180+1:end],t)
    write(io, "#Maxwell2\n")
    for i = 1:length(t)
        fx1 = Gx1[i]/Gmax
        fx2 = Gx2[i]/Gmax
        line = @sprintf "% .4f % .4f % .4f % .4f % .4f % .4f\n" fx1 fx1 fx1 fx2 fx2 fx2
        write(io, line)
    end
end;
## READ
f = readlines("./Maxwell2.txt")
G = zeros(size(f[2:end],1),6)
for (i, l) in enumerate(f[2:end])
    g = tryparse.(Float64, split(l," "))
    G[i,:] = g[g .≠ nothing]
end
using Plots
pyplot()
plot(G[:,1])
plot(G[:,4])