module MacroUnsafeBlink

"""
    @unsafe
A workaround for an upstream bug on Ubuntu 24. It disables sandboxing in Electron. 
Run this macro before calling `gogui()`.
The recommended way of using PackageMaker on Ubuntu 24 is however to use from VSCode, as 
in this case it works OK without calling `@unsafe`
# Examples
```julia-repl
julia> @unsafe;
julia> gogui()
``` 
"""
macro unsafe_blink()
    return quote
        Core.eval(@__MODULE__, quote
            function Blink.AtomShell.init(; debug = false)
                Blink.AtomShell.electron() # Check path exists
                p, dp = Blink.AtomShell.port(), Blink.AtomShell.port()
                debug && Blink.AtomShell.inspector(dp)
                dbg = debug ? "--debug=$dp" : []
                proc = (debug ? Blink.AtomShell.run_rdr : Blink.AtomShell.run)(
                    `$(Blink.AtomShell.electron()) --no-sandbox $dbg $(Blink.AtomShell.mainjs) port $p`; wait=false)
                conn = Blink.AtomShell.try_connect(ip"127.0.0.1", p)
                shell = Blink.AtomShell.Electron(proc, conn)
                Blink.AtomShell.initcbs(shell)
                return shell
            end
        end)
    end
end
export @unsafe_blink
end