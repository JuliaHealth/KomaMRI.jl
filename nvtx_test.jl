using Pkg
Pkg.activate(".")

using KomaMRI  

phantom = brain_phantom2D()[1:100]

if ARGS[1] == "simple"
    # Simple Motion
    print("Simple Motion\n")
    phantom.motion = SimpleMotion(
        ux = (x,y,z,t)->sin.(t),
        uy = (x,y,z,t)->sin.(t),
        uz = (x,y,z,t)->0
    )
elseif ARGS[1] == "arbitrary"
    # ArbitraryMotion
    print("Arbitrary Motion\n")
    K = 10
    Ns = length(phantom.x)
    phantom.motion = ArbitraryMotion(
        [1.0],
        K,
        rand(Ns, K-1),
        rand(Ns, K-1),
        zeros(Ns, K-1),
        Bool.(zeros(Ns, K-1))
    )
end

global sys = Scanner()

global seq = PulseDesigner.GRE(0.04,    # FOV
                               64,      # N
                               5e-3,    # TE
                               20e-3,    # TR
                               15,       # α
                               sys;
                               Δf = 0)

## Simulation
simulate(phantom, seq, sys)