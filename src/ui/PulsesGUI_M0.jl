plt = Observable{Any}(plot_M0(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_M0(p;darkmode), plt, seq_obs)
@js_ w document.getElementById("content").dataset.content = $CONT_M0
content!(w, "div#content", ui)
