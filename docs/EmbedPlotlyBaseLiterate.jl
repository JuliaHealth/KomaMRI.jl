import Base64: base64encode
import PlotlyBase

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

# This is type piracy, but it only affects the docs.
function Base.show(io::IO, ::MIME"text/html", fig::PlotlyBase.Plot)
    # Copy is required because we pop layout width/height for Plotly's inner JSON.
    # Mutating fig directly would change its user-visible behavior.
    plot = copy(fig)
    plot.frames = fig.frames
    plot.config = deepcopy(fig.config)
    plot.config.displayModeBar = false
    default_width = layout_to_html_default_size!(plot.layout.fields, :width, get(plot.layout.fields, :width, nothing))
    default_height = layout_to_html_default_size!(plot.layout.fields, :height, get(plot.layout.fields, :height, nothing))

    html_buffer = IOBuffer()
    PlotlyBase.to_html(
        html_buffer,
        plot;
        autoplay=false,
        full_html=true,
        include_plotlyjs="cdn",
        default_width=default_width,
        default_height=default_height,
    )

    html = String(take!(html_buffer))
    html = replace(html, "<body>" => "<body style=\"margin:0;overflow:hidden;background:transparent;\">")
    meta = get(plot.layout, :meta, nothing)
    koma = isnothing(meta) ? nothing : get(meta, :koma, nothing)
    if !isnothing(koma) && get(koma, :loop, false)
        html = replace(
            html,
            "</body>" => """
            <script>
            (() => {
                const plot = document.getElementById('$(plot.divid)');
                const connect = () => {
                    if (typeof plot.on !== 'function') return requestAnimationFrame(connect);
                    const play = plot.layout.updatemenus[0].buttons[0];
                    let looping = false;
                    plot.on('plotly_buttonclicked', event =>
                        looping = event.button.label === 'Play'
                    );
                    plot.on('plotly_animated', () => looping && requestAnimationFrame(() =>
                        looping && Plotly.animate(plot, null, play.args[1])
                    ));
                };
                connect();
            })();
            </script>
            </body>
            """,
        )
    end
    encoded = base64encode(html)
    onload = replace(
        """
        this.onload=null;
        this.removeAttribute('onload');
        const encoded=this.dataset.plotlyHtml;
        this.removeAttribute('data-plotly-html');
        const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));
        const doc=this.contentDocument;
        doc.open();
        doc.write(new TextDecoder().decode(bytes));
        doc.close();
        """,
        '&' => "&amp;",
        '"' => "&quot;",
        '<' => "&lt;",
        '>' => "&gt;",
    )
    width = to_css_size(get(fig.layout.fields, :width, nothing))
    _, plot_height = size(fig)
    height = to_css_size(plot_height)
    write(
        io,
        """
        <object type="text/html" data="about:blank" data-plotly-html="$encoded" onload="$onload" style="display:block;max-width:none;border:0;width:$width;height:$height;"></object>
        """,
    )
end
