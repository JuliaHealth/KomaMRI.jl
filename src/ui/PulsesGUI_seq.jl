#plt = Observable{Any}(plot_seq(seq; darkmode, range=[0 30]))
#ui = dom"div"(plt)
#map!(p->plot_seq(p; darkmode, range=[0 30]), plt, seq_obs)
#content!(w, "div#content", ui)

plt = Observable{Any}(plot_seq(seq; darkmode, range=[0 30]))
btn = button("Export .mat")
ui = dom"div"(vbox(dom"div"(hbox(btn)), plt))
map!(p->plot_seq(p; darkmode, range=[0 30]), plt, seq_obs)
content!(w, "div#content", ui)

function export2mat()
    max_rf_samples=100
    N = length(seq)
    ΔT = KomaMRICore.durs(seq)
    T0 = cumsum([0; ΔT],dims=1)
    off_val = Inf #This removes the unnecessary points in the plot

    #GRADS
    t1x = vcat([KomaMRICore.get_theo_t(seq.GR[1,i]) .+ T0[i] for i=1:N]...)
    t1y = vcat([KomaMRICore.get_theo_t(seq.GR[2,i]) .+ T0[i] for i=1:N]...)
    t1z = vcat([KomaMRICore.get_theo_t(seq.GR[3,i]) .+ T0[i] for i=1:N]...)
    Gx =  vcat([KomaMRICore.get_theo_A(seq.GR[1,i];off_val) for i=1:N]...)
    Gy =  vcat([KomaMRICore.get_theo_A(seq.GR[2,i];off_val) for i=1:N]...)
    Gz =  vcat([KomaMRICore.get_theo_A(seq.GR[3,i];off_val) for i=1:N]...)
    GRADS = hcat(t1x, t1y, t1z, Gx, Gy, Gz)
    #RFS
    t2 =  vcat([KomaMRICore.get_theo_t(seq.RF[1,i];max_rf_samples) .+ T0[i] for i=1:N]...)
    R =   vcat([KomaMRICore.get_theo_A(r;off_val,max_rf_samples) for r = seq.RF]...)
    RFS = hcat(t2, R)
    #ADC
    t3 =  vcat([KomaMRICore.get_theo_t(seq.ADC[i])  .+ T0[i] for i=1:N]...)
    D =   vcat([KomaMRICore.get_theo_A(d;off_val) for d = seq.ADC]...)
    ADCS = hcat(t3, D)

    seq_dict = Dict("GRAD" => GRADS,
                    "RF" => RFS,
                    "ADC" => ADCS,
                    "DUR" => seq.DUR,
                    "DEF" => seq.DEF)

    matwrite("sequence.mat", Dict("sequence" => seq_dict))

end

on(n -> export2mat(), btn)
