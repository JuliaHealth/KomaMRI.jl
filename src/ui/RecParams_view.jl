ui = plot_dict(recParams)
title = """<h1 style="padding: 8px 16px; color: #868888;">Reconstruction parameters</h1>"""
@js_ w document.getElementById("content").dataset.content = "recparams"
content!(w, "div#content", title*ui)
