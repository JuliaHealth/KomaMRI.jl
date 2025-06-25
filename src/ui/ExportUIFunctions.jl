"""
Returns the blink window with some custom styles and js logic.
"""
function setup_blink_window(; darkmode=true, frame=true, dev_tools=false, show_window=true)
    komamri_src_ui = @__DIR__
    komamri_root = dirname(dirname(komamri_src_ui))
    # Asset folders
    assets = AssetRegistry.register(komamri_root * "/assets/")
    scripts = AssetRegistry.register(komamri_src_ui * "/scripts/")
    css = AssetRegistry.register(komamri_src_ui * "/css/")
    # Images
    logo = joinpath(assets, "logo-dark.svg")
    home_image = joinpath(assets, "home-image.svg")
    app_icon = komamri_root * "/assets/app-icon.png"
    # JS
    bsjs = joinpath(scripts, "bootstrap.bundle.min.js") #this already has Popper
    bscss = joinpath(css, "bootstrap.min.css")
    bsiconcss = joinpath(css, "bootstrap-icons.css")
    jquery = joinpath(scripts, "jquery-3.4.1.slim.min.js")
    # mathjaxsetup = joinpath(scripts, "mathjaxsetup.js")
    # KaTeX
    katexrender = joinpath(scripts, "auto-render.min.js")
    katexjs = joinpath(scripts, "katex.min.js")
    katexcss = joinpath(css, "katex.min.css")
    # User defined JS and CSS
    customcss = joinpath(css, "custom.css")
    customjs = joinpath(scripts, "custom.js")
    sidebarcss = joinpath(css, "sidebars.css")
    # Custom icons
    icons = joinpath(css, "icons.css")
    ## WINDOW
    w = Blink.Window(
        Dict(
            "title" => "KomaUI",
            "autoHideMenuBar" => true,
            "frame" => frame, #removes title bar
            "node-integration" => true,
            :icon => app_icon,
            "width" => 1200,
            "height" => 800,
            :show => show_window,
        );
        async=false,
    )

    ## NAV BAR
    sidebar = open(f -> read(f, String), komamri_src_ui * "/html/sidebar.html")
    sidebar = replace(sidebar, "LOGO" => logo)
    ## CONTENT
    index = open(f -> read(f, String), komamri_src_ui * "/html/index.html")
    index = replace(index, "ICON" => home_image)

    @sync begin
        ## CSS
        loadcss!(w, bscss)
        loadcss!(w, bsiconcss)
        loadcss!(w, customcss)
        loadcss!(w, sidebarcss)
        # KATEX
        loadcss!(w, katexcss)
        loadjs!(w, katexjs)
        loadjs!(w, katexrender)
        # JQUERY, BOOSTRAP JS
        loadjs!(w, customjs)    #must be before jquery
        loadjs!(w, jquery)
        loadjs!(w, bsjs)        #after jquery
        # LOAD ICONS
        loadcss!(w, icons)
    end

    #Update GUI's home
    body!(w, *(sidebar, index); async=false) #async false is important to set the background correctly
    if darkmode
        @js_ w document.getElementById("main").style = "background-color:rgb(13,16,17);"
    end
    # Return the Blink window
    if dev_tools
        Blink.tools(w)
    end
    return w, index
end

"""
Returns the default scanner used by the UI and print some information about it.
"""
function setup_scanner()

    # Print information and get the default Scanner struct
    @info "Loaded `Scanner` to `sys_ui[]`"
    sys = Scanner()

    # Return the default Scanner struct
    return sys
end

"""
Returns the default sequence used by the UI and print some information about it.
"""
function setup_sequence(sys::Scanner)

    # Print information and get the default Sequence struct
    @info "Loaded `Sequence` to `seq_ui[]`"
    seq = PulseDesigner.EPI_example(; sys)

    # Return the default Sequence struct
    return seq
end

"""
Returns the default phantom used by the UI and print some information about it.
"""
function setup_phantom(; phantom_mode="2D")

    # Print information and get the default Phantom struct
    @info "Loaded `Phantom` to `obj_ui[]`"
    obj = phantom_mode == "3D" ? brain_phantom3D() : brain_phantom2D()
    obj.Δw .= 0 # Removes off-resonance

    # Return the default Phantom struct
    return obj
end

"""
Returns the default raw signal used by the UI.
"""
function setup_raw()

    # Define the default RawAcquisitionData struct
    raw = RawAcquisitionData(
        Dict(
            "systemVendor" => "",
            "encodedSize" => [2, 2, 1],
            "reconSize" => [2, 2, 1],
            "number_of_samples" => 4,
            "encodedFOV" => [100.0, 100.0, 1],
            "trajectory" => "other",
        ),
        [
            KomaMRICore.Profile(
                AcquisitionHeader(; trajectory_dimensions=2, sample_time_us=1),
                [0.0 0.0 1 1; 0 1 1 1] ./ 2,
                reshape([0.0; 0im; 0; 0], 4, 1),
            ),
        ],
    )

    # Return the default RawAcquisitionData struct
    return raw
end

"""
Updates the blink window with a loading content
"""
function display_loading!(w::Window, msg::String)
    komamri_src_ui = @__DIR__
    loading = replace(
        open((f) -> read(f, String), komamri_src_ui * "/html/loading.html"),
        "LOADDES" => msg,
    )
    return content!(w, "div#content", loading)
end

"""
Fileficker callback when sequence is loaded
"""
function callback_filepicker(filename::String, w::Window, seq::Sequence)
    if filename != ""
        filename_extension = splitext(filename)[end]
        if filename_extension == ".seqk"     # Koma
            seq = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(filename), "seq")
        elseif filename_extension == ".seq" # Pulseq
            seq = read_seq(filename)        # Pulseq read
        end
        @js_ w (
            @var name = $(basename(filename));
            document.getElementById("seqname").innerHTML =
                "<abbr title='" + name + "'>" + name + "</abbr>";
            Toasty(
                "0",
                "Loaded <b>" + name + "</b> successfully",
                """
<ul>
    <li>
        <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('pulses_seq', 1)"><i class="fa fa-search"></i></button>
        Updating <b>Sequence</b> plots ...
    </li>
    <li>
        <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('simulate', 1)"><i class="bi bi-caret-right-fill"></i></button>
        Ready to <b>simulate</b>?
    </li>
</ul>
""",
            )
        )
        return seq
    end
    return seq
end

"""
Fileficker callback when phantom is loaded
"""
function callback_filepicker(filename::String, w::Window, obj::Phantom)
    if filename != ""
        filename_extension = splitext(filename)[end]
        if filename_extension == ".phantom" # Koma
            obj = read_phantom(filename)
        elseif filename_extension == ".h5"  # JEMRIS
            obj = read_phantom_jemris(filename)
        end
        @js_ w (
            @var name = $(basename(filename));
            document.getElementById("phaname").innerHTML =
                "<abbr title='" + name + "'>" + name + "</abbr>";
            Toasty(
                "0",
                "Loaded <b>" + name + "</b> successfully",
                """
<ul>
    <li>
        <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('phantom', 1)"><i class="fa fa-search"></i></button>
        Updating <b>Phantom</b> plots ...
    </li>
    <li>
        <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('simulate', 1)"><i class="bi bi-caret-right-fill"></i></button>
        Ready to <b>simulate</b>?
    </li>
</ul>
""",
            )
        )
        return obj
    end
    return obj
end

"""
Fileficker callback when raw signal is loaded
"""
function callback_filepicker(filename::String, w::Window, raw::RawAcquisitionData)
    if filename != ""
        raw = RawAcquisitionData(ISMRMRDFile(filename))
        if raw.params["systemVendor"] != "KomaMRI.jl"
            @warn "ISMRMRD files generated externally could cause problems during the reconstruction. We are currently improving compatibility."
        end
        @js_ w (
            @var name = $(basename(filename));
            document.getElementById("rawname").innerHTML =
                "<abbr title='" + name + "'>" + name + "</abbr>";
            Toasty(
                "0",
                "Loaded <b>" + name + "</b> successfully",
                """
<ul>
    <li>
        <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('sig', 1)"><i class="fa fa-search"></i></button>
        Updating <b>Raw data</b> plots ...
    </li>
    <li>
        <button class="btn btn-success btn-circle btn-circle-sm m-1" onclick="Blink.msg('recon', 1)"><i class="bi bi-caret-right-fill"></i></button>
        Ready to <b>reconstruct</b>?
    </li>
</ul>
""",
            )
        )
        return raw
    end
    return raw
end

"""
Updates the UI with sequence plots
"""
function view_ui!(seq::Sequence, w::Window; type="sequence", darkmode=true)
    # Add loading progress and then a plot to the UI depending on type of the plot
    if type == "sequence"
        display_loading!(w, "Plotting sequence ...")
        long_seq = length(seq) > 1_000
        time_end = !long_seq ? 30 : dur(seq) * 1e3
        widget_plot = plot_seq(seq; darkmode, range=[0 time_end], slider=!long_seq, gl=long_seq, show_adc=false)
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "sequence"
    elseif type == "kspace"
        display_loading!(w, "Plotting kspace ...")
        widget_plot = plot_kspace(seq; darkmode)
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "kspace"
    elseif type == "moment0"
        display_loading!(w, "Plotting moment 0 ...")
        widget_plot = plot_M0(seq; darkmode)
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "m0"
    elseif type == "moment1"
        display_loading!(w, "Plotting moment 1 ...")
        widget_plot = plot_M1(seq; darkmode)
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "m1"
    elseif type == "moment2"
        display_loading!(w, "Plotting moment 2 ...")
        widget_plot = plot_M2(seq; darkmode)
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "m2"
    end
end

"""
Updates the UI with phantom plots
"""
function view_ui_phantom!(
    obj::Phantom,
    w::Window,
    seq::Sequence,
    buttons_obj::Vector{Widget{:button,Int64}};
    key=:ρ,
    darkmode=true,
)
    display_loading!(w, "Plotting phantom ...")
    widget_plot = plot_phantom_map(obj, key; time_samples=5, darkmode)
    div_content = dom"div"(hbox(buttons_obj...), widget_plot)
    content!(w, "div#content", div_content)
    @js_ w document.getElementById("content").dataset.content = "phantom"
end
function view_ui!(
    obj::Phantom,
    w::Window,
    seq::Sequence,
    buttons_obj::Vector{Widget{:button,Int64}};
    key=:ρ,
    darkmode=true,
)
    return view_ui_phantom!(obj, w, seq, buttons_obj; key, darkmode)
end
function view_ui!(
    cnt::Integer,
    w::Window,
    obj::Phantom,
    seq::Sequence,
    buttons_obj::Vector{Widget{:button,Int64}};
    key=:ρ,
    darkmode=true,
)
    return view_ui_phantom!(obj, w, seq, buttons_obj; key, darkmode)
end

"""
Updates the UI with scanner information
"""
function view_ui!(sys::Scanner, w::Window)
    display_loading!(w, "Displaying scanner parameters ...")
    sys_dict = Dict(
        "B0" => sys.B0,
        "B1" => sys.B1,
        "Gmax" => sys.Gmax,
        "Smax" => sys.Smax,
        "ADC_dt" => sys.ADC_Δt,
        "seq_dt" => sys.seq_Δt,
        "GR_dt" => sys.GR_Δt,
        "RF_dt" => sys.RF_Δt,
        "RF_ring_down_T" => sys.RF_ring_down_T,
        "RF_dead_time_T" => sys.RF_dead_time_T,
        "ADC_dead_time_T" => sys.ADC_dead_time_T,
    )
    plt = plot_dict(sys_dict)
    title = """<h1 style="padding: 8px 16px; color: #868888;">Scanner parameters</h1>"""
    content!(w, "div#content", title * plt)
    @js_ w document.getElementById("content").dataset.content = "scanneparams"
end

"""
Updates the UI with simulation parameters information
"""
function view_ui!(sim_params::Dict{String,Any}, w::Window)
    display_loading!(w, "Displaying simulation parameters ...")
    plt = plot_dict(sim_params)
    title = """<h1 style="padding: 8px 16px; color: #868888;">Simulation parameters</h1>"""
    content!(w, "div#content", title * plt)
    @js_ w document.getElementById("content").dataset.content = "simparams"
end

"""
Updates the UI with raw signal plot
"""
function view_ui!(raw::RawAcquisitionData, w::Window; darkmode=true)
    display_loading!(w, "Plotting raw signal ...")
    widget_plot = plot_signal(raw; darkmode)
    content!(w, "div#content", dom"div"(widget_plot))
    @js_ w document.getElementById("content").dataset.content = "sig"
end

"""
Updates the UI with reconstruction parameters information
"""
function view_ui!(rec_params::Dict{Symbol,Any}, w::Window)
    display_loading!(w, "Displaying reconstruction parameters ...")
    plt = plot_dict(rec_params)
    title = """<h1 style="padding: 8px 16px; color: #868888;">Reconstruction parameters</h1>"""
    content!(w, "div#content", title * plt)
    @js_ w document.getElementById("content").dataset.content = "recparams"
end

"""
Updates the UI with image plots
"""
function view_ui!(img::Array, w::Window; type="absi", darkmode=true)
    # Add loading progress and then a plot to the UI depending on type of the plot
    if type == "absi"
        display_loading!(w, "Plotting image magnitude ...")
        widget_plot = @manipulate for slice in 1:size(img, 3)
            aux = abs.(img) * prod(size(img)[1:2])
            plot_image(
                aux[:, :, slice];
                zmin=minimum(aux[:]),
                zmax=maximum(aux[:]),
                darkmode,
                title="Reconstruction ($slice/$(size(img, 3)))",
            )
        end
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "absi"
    elseif type == "angi"
        display_loading!(w, "Plotting image phase ...")
        widget_plot = @manipulate for slice in 1:size(img, 3)
            aux = angle.(img[:, :, slice])
            plot_image(
                aux;
                zmin=-π,
                zmax=π,
                darkmode,
                title="Reconstruction ($slice/$(size(img, 3)))",
            )
        end
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "angi"
    elseif type == "absk"
        display_loading!(w, "Plotting image k ...")
        widget_plot = @manipulate for slice in 1:size(img, 3)
            kspace = fftc(img)
            aux = log.(abs.(kspace[:, :, slice]) .+ 1)
            plot_image(
                aux;
                zmin=0,
                zmax=0.1 * maximum(aux[:]),
                darkmode,
                title="Reconstruction ($slice/$(size(img, 3)))",
            )
        end
        content!(w, "div#content", dom"div"(widget_plot))
        @js_ w document.getElementById("content").dataset.content = "absk"
    end
end

"""
Saves to .mat in the UI (displays a Toast message en the UI)
"""
function save_ui!(
    w::Window,
    seq::Sequence,
    obj::Phantom,
    sys::Scanner,
    raw,
    img,
    rec_params,
    mat_folder;
    type="all",
)
    if type == "all"
        str_toast = export_2_mat(
            seq, obj, sys, raw, rec_params, img, mat_folder; type, matfilename=""
        )
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfolder"
    elseif type == "sequence"
        str_toast = export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfolderseq"
    elseif type == "phantom"
        str_toast = export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfolderpha"
    elseif type == "scanner"
        str_toast = export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfoldersca"
    elseif type == "raw"
        str_toast = export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfolderraw"
    elseif type == "image"
        str_toast = export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
        @js_ w (@var msg = $str_toast; Toasty("1", "Saved .mat files", msg))
        @js_ w document.getElementById("content").dataset.content = "matfolderima"
    end
end
