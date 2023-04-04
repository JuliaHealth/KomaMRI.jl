#plt = Observable{Any}(plot_kspace(seq;darkmode))
#ui = dom"div"(plt)
#map!(p->plot_kspace(p;darkmode), plt, seq_obs)
#content!(w, "div#content", ui)

plt = Observable{Any}(plot_kspace(seq;darkmode))
btn = button("Export .mat")
#ui = dom"div"(vbox(dom"div"(hbox(btn)), plt))
ui = dom"div"(plt)
map!(p->plot_kspace(p;darkmode), plt, seq_obs)
content!(w, "div#content", ui)

function export2mat()

    kspace, kspace_adc = get_kspace(seq; Δt=1)

    matwrite("kspace.mat", Dict("kspace" => kspace, "kspace_adc" => kspace_adc))

end

on(n -> export2mat(), btn)
