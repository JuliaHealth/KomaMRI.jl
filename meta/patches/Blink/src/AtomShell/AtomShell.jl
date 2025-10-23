module AtomShell

using ..Blink: resource, @init, @errs, get_blink_root_dir
using Sockets
using WebIO

abstract type Shell end

include("process.jl")
include("window.jl")
include("webio.jl")

end
