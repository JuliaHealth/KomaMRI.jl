plt = Observable{Any}(plot_M2(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_M2(p;darkmode), plt, seq_obs)
@js_ w document.getElementById("content").dataset.content = $CONT_M2
content!(w, "div#content", ui)
