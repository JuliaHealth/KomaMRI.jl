plt = Observable{Any}(plot_seq(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_seq(p;darkmode), plt, seq_obs)
content!(w, "div#content", ui)