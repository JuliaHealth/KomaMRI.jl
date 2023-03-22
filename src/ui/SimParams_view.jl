#ui = plot_dict(simParams)
#title = """<h1 class="text-white">Simulation parameters</h1>"""
#content!(w, "div#content", title*ui)

plt = plot_dict(simParams)
title = """<h1 class="text-white">Simulation parameters</h1>"""
btn = """<button class="btn btn-dark btn-sm" onclick='Blink.msg("btnsimparams", 1)'>Export .mat</button>"""
content!(w, "div#content", btn*title*plt)

handle(w, "btnsimparams") do args...
    matwrite("sim_params.mat", Dict("sim_params" => raw_ismrmrd.params["userParameters"]))
end
