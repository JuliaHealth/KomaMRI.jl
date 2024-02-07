using Pkg
Pkg.activate(".")

using KomaMRI  

KomaMRICore.CUDA.allowscalar(false)

phantom = brain_phantom2D()[1:100]

if ARGS[1] == "simple"
    # Simple Motion
    print("Simple Motion\n")
    phantom.motion = SimpleMotion(Translation(direction = [1.0,0,0], v = 0.1))

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

sys = Scanner()

seq = PulseDesigner.EPI_example(; sys)

## Simulation
simulate(phantom, seq, sys)
simulate(phantom, seq, sys)