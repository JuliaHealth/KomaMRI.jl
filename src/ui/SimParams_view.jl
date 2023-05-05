plt = plot_dict(simParams)
title = """<h1 style="padding: 8px 16px; color: #868888;">Simulation parameters</h1>"""
content!(w, "div#content", title*plt)
