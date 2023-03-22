#plt = Observable{Any}(plot_M0(seq;darkmode))
#ui = dom"div"(plt)
#map!(p->plot_M0(p;darkmode), plt, seq_obs)
#content!(w, "div#content", ui)

plt = Observable{Any}(plot_M0(seq;darkmode))
btn = button("Export .mat")
ui = dom"div"(vbox(dom"div"(hbox(btn)), plt))
map!(p->plot_M0(p;darkmode), plt, seq_obs)
content!(w, "div#content", ui)

function export2mat()

    dt = 1
    t, Δt = KomaMRICore.get_uniform_times(seq, dt)
    ts = t .+ Δt
    k, _ =  KomaMRICore.get_kspace(seq; Δt=dt)
    momentum0 = hcat(t, k)

    matwrite("momentum0.mat", Dict("momentum0" => momentum0))

end

on(n -> export2mat(), btn)
