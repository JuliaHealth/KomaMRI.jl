import Base64: base64encode
import PlotlyJS

# Convert Plotly width/height values to valid CSS size values.
to_css_size(x::Integer) = string(x, "px")
to_css_size(x::AbstractFloat) = string(round(Int, x), "px")
to_css_size(x::AbstractString) = x
to_css_size(::Nothing) = "100%"
to_css_size(x) = throw(ArgumentError("Unsupported Plotly size type: $(typeof(x))"))

# For CSS string sizes, keep the explicit size on the outer <object> only.
# The inner Plotly HTML should stay at 100% to avoid double-scaling.
layout_to_html_default_size!(fields, key::Symbol, value::AbstractString) = (pop!(fields, key, nothing); "100%")
layout_to_html_default_size!(fields, key::Symbol, ::Nothing) = (pop!(fields, key, nothing); "100%")
layout_to_html_default_size!(fields, key::Symbol, value::Number) = to_css_size(value)
layout_to_html_default_size!(fields, key::Symbol, value) = throw(ArgumentError("Unsupported Plotly layout $key type: $(typeof(value))"))

# This is type piracy! but it only affects the docs.
# We are extending the `Base.show` method for the `PlotlyJS.SyncPlot` 
function Base.show(io::IO, ::MIME"text/html", fig::PlotlyJS.SyncPlot)
    # Copy is required because we pop layout width/height for Plotly's inner JSON.
    # Mutating fig.plot directly would change size(fig) and user-visible behavior.
    plot = copy(fig.plot)
    default_width = layout_to_html_default_size!(plot.layout.fields, :width, get(plot.layout.fields, :width, nothing))
    default_height = layout_to_html_default_size!(plot.layout.fields, :height, get(plot.layout.fields, :height, nothing))

    html_buffer = IOBuffer()
    PlotlyJS.PlotlyBase.to_html(
        html_buffer,
        plot;
        full_html=true,
        include_plotlyjs="cdn",
        default_width=default_width,
        default_height=default_height,
    )

    html = String(take!(html_buffer))
    html = replace(html, "<body>" => "<body style=\"margin:0;overflow:hidden;\">")
    encoded = base64encode(html)
    width = to_css_size(get(fig.plot.layout.fields, :width, nothing))
    _, plot_height = size(fig)
    height = to_css_size(plot_height)
    write(
        io,
        """
        <object type="text/html" data="data:text/html;base64,$encoded" style="display:block;max-width:none;border:0;width:$width;height:$height;"></object>
        """,
    )
end
