module PlotlyKaleido

using JSON: JSON
using Artifacts: @artifact_str
using Base64: Base64
using Kaleido_jll: Kaleido_jll

export savefig

#-----------------------------------------------------------------------------# Windows Fallback

FALLBACK_DIR() = @static Sys.iswindows() ? artifact"Kaleido_fallback" : ""

get_kaleido_version() = read(joinpath(Kaleido_jll.artifact_dir, "version"), String)
should_try_fallback() = Sys.iswindows() && (get_kaleido_version() !== "0.1.0")

const USE_KALEIDO_FALLBACK = Ref(should_try_fallback())

#-----------------------------------------------------------------------------# Kaleido Process
mutable struct Pipes
    stdin::Pipe
    stdout::Pipe
    stderr::Pipe
    proc::Base.Process
    Pipes() = new()
end

const P = Pipes()

const _mathjax_url_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax"
const _mathjax_last_version = v"2.7.9"

function warn_and_kill(s::String)
    @warn "$s"
    kill_kaleido()
    return nothing
end

kill_kaleido() = is_running() && (kill(P.proc); wait(P.proc))

is_running() = isdefined(P, :proc) && isopen(P.stdin) && process_running(P.proc)

restart(; kwargs...) = (kill_kaleido(); start(; kwargs...))

# The content of this function is inspired from https://discourse.julialang.org/t/readline-with-default-value-if-no-input-after-timeout/100388/2?u=disberd
function readline_noblock(io; timeout = 10)
    msg = Channel{String}(1)

    task = Task() do
        try
            put!(msg, readline(io))
        catch
            put!(msg, "Stopped")
        end
    end

    interrupter = Task() do
        sleep(timeout)
        if !istaskdone(task)
            Base.throwto(task, InterruptException())
        end
    end

    schedule(interrupter)
    schedule(task)
    wait(task)
    out = take!(msg)
    if out === "Stopped" 
        warn_str = "It looks like the Kaleido process is not responding since $(timeout) seconds. 
The unresponsive process will be killed, but this means that you will not be able to save figures using `savefig`."

        if should_try_fallback() && !USE_KALEIDO_FALLBACK[]
            warn_str *= "

You seem to be on Windows but have disabled the automatic fallback to version 0.1 of Kaleido. You may want to try enabling it by calling `PlotlyKaleido.USE_KALEIDO_FALLBACK[] = true`, as higher version of Kaleido are known to have issues on Windows.
Check the Package Readme at https://github.com/JuliaPlots/PlotlyKaleido.jl/tree/main#windows-note for more details."
        end

        warn_str *= "

Alternatively, you might try using a longer timeout to check if the process is not responding by passing the desired value in seconds using the `timeout` kwarg when calling `PlotlyKaleido.start` or `PlotlyKaleido.restart`"
        warn_and_kill(warn_str)
    end
    return out
end

function get_base_cmd()
    cmd = if should_try_fallback() && USE_KALEIDO_FALLBACK[]
        # For the fallback we don't fully reproduce the jll machinery as this is much simpler and should work fine for kaleido specifically on windows.
        dir = FALLBACK_DIR()

        Cmd(`$(joinpath(dir, "bin", "kaleido.exe"))`; dir)
    else
        dir = Kaleido_jll.artifact_dir
        Cmd(Kaleido_jll.kaleido(); dir)
    end
    return cmd
end

function start(;
    plotly_version = missing,
    mathjax = missing,
    mathjax_version::VersionNumber = _mathjax_last_version,
    timeout = 10,
    kwargs...,
)
    is_running() && return
    # The kaleido executable must be run from the artifact directory
    BIN = get_base_cmd()
    # We push the mandatory plotly flag
    push!(BIN.exec, "plotly")
    chromium_flags = ["--disable-gpu", Sys.isapple() ? "--single-process" : "--no-sandbox"]
    extra_flags = if plotly_version === missing
        (; plotlyjs = string("file://", artifact"plotly/plotly.js"))
    else
        # We create a plotlyjs flag pointing at the specified plotly version
        (; plotlyjs = "https://cdn.plot.ly/plotly-$(plotly_version).min.js", kwargs...)
    end
    if mathjax === missing
        mathjaxfile = string("file://", artifact"mathjax/mathjax.js")
        push!(chromium_flags, "--mathjax=$mathjaxfile")
    else
        if mathjax_version > _mathjax_last_version
            error(
                "The given mathjax version ($(mathjax_version)) is greater than the last supported version ($(_mathjax_last_version)) of Kaleido.",
            )
        end
        if mathjax isa Bool
            mathjax && push!(
                chromium_flags,
                "--mathjax=$(_mathjax_url_path)/$(mathjax_version)/MathJax.js",
            )
        elseif mathjax isa String
            # We expect the keyword argument to be a valid URL or similar, else error "Kaleido startup failed with code 1".
            push!(chromium_flags, "--mathjax=$(mathjax)")
        else
            @warn """The value of the provided argument
                    mathjax=$(mathjax)
                  is neither a Bool nor a String and has been ignored."""
        end
    end
    # Taken inspiration from https://github.com/plotly/Kaleido/blob/3b590b563385567f257db8ff27adae1adf77821f/repos/kaleido/py/kaleido/scopes/base.py#L116-L141
    user_flags = String[]
    for (k, v) in pairs(extra_flags)
        flag_name = replace(string(k), "_" => "-")
        if v isa Bool
            v && push!(user_flags, "--$flag_name")
        else
            push!(user_flags, "--$flag_name=$v")
        end
    end
    # We add the flags to the BIN
    append!(BIN.exec, chromium_flags, user_flags)

    kstdin = Pipe()
    kstdout = Pipe()
    kstderr = Pipe()
    kproc = 
        run(pipeline(BIN, stdin = kstdin, stdout = kstdout, stderr = kstderr), wait = false)

    process_running(kproc) || error("There was a problem starting up kaleido.")
    close(kstdout.in)
    close(kstderr.in)
    close(kstdin.out)
    Base.start_reading(kstderr.out)

    global P
    P.stdin = kstdin
    P.stdout = kstdout
    P.stderr = kstderr
    P.proc = kproc

    res = readline_noblock(P.stdout; timeout)  # {"code": 0, "message": "Success", "result": null, "version": "0.2.1"}
    length(res) == 0 && warn_and_kill("Kaleido startup failed.")
    if is_running()
        code = JSON.parse(res)["code"]
        code == 0 || warn_and_kill("Kaleido startup failed with code $code.")
    end
    return
end


#-----------------------------------------------------------------------------# save
const ALL_FORMATS = ["png", "jpeg", "webp", "svg", "pdf", "eps", "json"]
const TEXT_FORMATS = ["svg", "json", "eps"]


function save_payload(io::IO, payload::AbstractString, format::AbstractString)
    is_running() || error("It looks like the Kaleido process is not running, so you can not save plotly figures.
Remember to start the process before using `savefig` by calling `PlotlyKaleido.start()`.
If the process was killed due to an error during initialization, you will receive a warning when the `PlotlyKaleido.start` function is executing")
    format in ALL_FORMATS || error("Unknown format $format. Expected one of $ALL_FORMATS")

    bytes = transcode(UInt8, payload)
    write(P.stdin, bytes)
    write(P.stdin, transcode(UInt8, "\n"))
    flush(P.stdin)

    res = readline(P.stdout)
    obj = JSON.parse(res)
    obj["code"] == 0 || error("Transform failed: $res")

    img = String(obj["result"])

    # base64 decode if needed, otherwise transcode to vector of byte
    bytes = format in TEXT_FORMATS ? transcode(UInt8, img) : Base64.base64decode(img)

    write(io, bytes)
end

function savefig(io::IO, plot; height = 500, width = 700, scale = 1, format = "png")
    payload = JSON.json((; height, width, scale, format, data = plot))
    save_payload(io, payload, format)
end

function savefig(
    io::IO,
    plot::AbstractString;
    height = 500,
    width = 700,
    scale = 1,
    format = "png",
)
    payload = "{\"width\":$width,\"height\":$height,\"scale\":$scale,\"data\": $plot}"
    save_payload(io, payload, format)
end

function savefig(filename::AbstractString, plot; kw...)
    format = get(kw, :format, split(filename, '.')[end])
    open(io -> savefig(io, plot; format, kw...), filename, "w")
    filename
end

savefig(plot, filename::AbstractString; kw...) = savefig(filename, plot; kw...)

end # module Kaleido
