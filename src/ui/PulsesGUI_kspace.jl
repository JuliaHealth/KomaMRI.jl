plt = Observable{Any}(plot_kspace(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_kspace(p;darkmode), plt, seq_obs)
content!(w, "div#content", ui)
