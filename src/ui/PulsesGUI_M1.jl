plt = Observable{Any}(plot_M1(seq;darkmode))
ui = dom"div"(plt)
map!(p->plot_M1(p;darkmode), plt, seq_obs)
@js_ w document.getElementById("content").dataset.content = $CONT_M1
content!(w, "div#content", ui)
