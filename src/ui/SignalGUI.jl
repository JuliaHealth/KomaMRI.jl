plt = Observable{Any}(plot_signal(raw_ismrmrd; darkmode))
ui = dom"div"(plt)
map!(p->plot_signal(p; darkmode), plt, sig_obs)
content!(w, "div#content", ui)
@js_ w document.getElementById("content").dataset.content = "sig"
