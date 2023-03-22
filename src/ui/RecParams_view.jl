ui = plot_dict(recParams)
title = """<h1 class="text-white">Reconstruction parameters</h1>"""
btn = """<button class="btn btn-dark btn-sm" onclick='Blink.msg("btnrecparams", 1)'>Export .mat</button>"""
content!(w, "div#content", btn*title*ui)

handle(w, "btnrecparams") do args...
    recParams_dict = Dict("reco" => recParams[:reco],
                        "Nx" => recParams[:reconSize][1],
                        "Ny" => recParams[:reconSize][2])
    matwrite("rec_params.mat", Dict("rec_params" => recParams_dict))
end
