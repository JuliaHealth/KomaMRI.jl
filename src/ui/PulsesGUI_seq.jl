plt = Observable{Any}(plot_seq(seq; darkmode, range=[0 20]))
ui = dom"div"(plt)
map!(p->plot_seq(p; darkmode, range=[0 20]), plt, seq_obs)
content!(w, "div#content", ui)