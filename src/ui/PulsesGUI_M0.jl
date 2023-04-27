plt = Observable{Any}(plot_M0(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_M0(p;darkmode), plt, seq_obs)
content!(w, "div#content", ui)
