__precompile__()

module Blink

using Reexport
using Distributed: Future
using Sockets
using Logging
using Base64: stringmime
using WebIO

function get_blink_root_dir()
    
    dir = Base.pkgdir(@__MODULE__)

    if isnothing(dir)
        # An ugly hack
        webio_dir = Base.pkgdir(WebIO)
        
        if isnothing(webio_dir)
            return @show joinpath(@__DIR__, "..")
        else
            blink_dir = joinpath(dirname(webio_dir), "Blink")
            return @show blink_dir
        end
    else
        return @show dir
    end
end

include("lazy/lazy.jl")

@init begin
    @show get_blink_root_dir()
    @show WebIO.get_webio_root_dir()
    WebIO.__init__()
    
end

include("rpc/rpc.jl")
include("content/content.jl")

include("AtomShell/AtomShell.jl")
export AtomShell
@reexport using .AtomShell
import .AtomShell: resolve_blink_asset



end # module
