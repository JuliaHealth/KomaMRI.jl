using MRIsim, JuMP, Ipopt
using MRIsim: plot_grads, get_Bmatrix, get_SRmatrix, get_Mmatrix, 
              Grad_fun, dur, plot_grads_moments, get_bvalue
#Parameters
Gmax =  30e-3 #21.7e-3; 
Smax = 52.5
τ = 100e-3
Nt = 600
DIF = Sequence([Grad_fun(x-> 1, τ, Nt); Grad_fun(x-> 1, τ, Nt)])
_, N = size(DIF.GR)
B, SR, M =  get_Bmatrix(DIF)[1,:,:], get_SRmatrix(DIF), get_Mmatrix(DIF)
start180, end180 = 400-50, 400+50 
DIF = [copy(DIF) for k=0:3]
for k = 0:3
    Mm = M[1:k+1,:]
    #optimization
    model = Model(); set_optimizer(model, Ipopt.Optimizer); set_silent(model)
    @variable(model, -Gmax <= x[1:N] <= Gmax, start=1); #max-grads
    @objective(model, Min, -x'*B*x); #b-value
    @constraint(model, rf180,  x[start180:end180] .== 0); #rf
    @constraint(model, slewrate, -Smax .<= SR*x .<= Smax); #slewrate
    @constraint(model, moments, Mm*x .== 0); #moments
    # @constraint(model, concomitant, sum(x[1:start180-1].^2) == sum(x[end180+1:end].^2)); #concomitant
    optimize!(model);
    gx = value.(x) #retrieving solution
    #seq object
    setproperty!.(DIF[k+1].GR, :A, [gx'; 0*gx']);
    @info "$k-order motion-compensated"
    @info "b-value: $(round(gx'*B*gx, digits=2)) s/mm2"
end
##Plots
p1 = plot_grads_moments(DIF[1],title="Maxwell,M0-compensated, b=$(round(get_bvalue(DIF[1])*1e-6, digits=2)) s/mm2")
p2 = plot_grads_moments(DIF[2],title="Maxwell,M0,M1-compensated, b=$(round(get_bvalue(DIF[2])*1e-6, digits=2)) s/mm2")
p3 = plot_grads_moments(DIF[3],title="Maxwell,M0,M1,M2-compensated, b=$(round(get_bvalue(DIF[3])*1e-6, digits=2)) s/mm2")
p4 = plot_grads_moments(DIF[4],title="Maxwell,M0,M1,M2,M3-compensated, b=$(round(get_bvalue(DIF[4])*1e-6, digits=2)) s/mm2")
##
using PlotlyJS
PlotlyJS.savefig(p1, "./src/examples/Figures/M0.pdf", width=1000)
PlotlyJS.savefig(p2, "./src/examples/Figures/M1.pdf", width=1000)
PlotlyJS.savefig(p3, "./src/examples/Figures/M2.pdf", width=1000)
PlotlyJS.savefig(p4, "./src/examples/Figures/M3.pdf", width=1000)
