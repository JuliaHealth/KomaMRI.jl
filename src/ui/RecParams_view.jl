ui = plot_dict(recParams)
title = """<h1 style="padding: 8px 16px;">Reconstruction parameters</h1>"""
content!(w, "div#content", title*ui)
