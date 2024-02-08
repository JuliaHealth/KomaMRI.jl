using Pkg
Pkg.activate(".")

using KomaMRI  

# KomaMRICore.CUDA.allowscalar(true)

phantom = brain_phantom2D()
phantom = phantom[1:minimum([length(phantom.x),parse(Int64,ARGS[1])])]
print("Spins: ", length(phantom.x), "\n")

if ARGS[2] == "simple"
    # Simple Motion
    print("Simple Motion\n")
    phantom.motion = SimpleMotion(Translation([1.0,0,0],0.1))
    # phantom.motion = SimpleMotion(Rotation([0,1,0],[0,0,0],1))

elseif ARGS[2] == "arbitrary"
    # ArbitraryMotion
    print("Arbitrary Motion\n")
    K = 2
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
sleep(1)
simulate(phantom, seq, sys)