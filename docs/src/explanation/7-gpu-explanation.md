# GPU Parallelization

KomaMRI uses a vendor agnostic approach to GPU parallelization in order to support multiple GPU backends. Currently, the following backends are supported:

* CUDA.jl (Nvidia)
* Metal.jl (Apple)
* AMDGPU.jl (AMD)
* oneAPI.jl (Intel)

## Choosing a GPU Backend

To determine which backend to use, KomaMRI uses [package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) (introduced in Julia 1.9) to avoid having the packages for each GPU backend as explicit dependencies. This means that the user is responsible for loading the backend package (e.g. `using CUDA`) at the beginning of their code, or prior to calling KomaUI(), otherwise, Koma will default back to the CPU:

```julia
using KomaMRI
using CUDA # loading CUDA will load KomaMRICoreCUDAExt, selecting the backend
```

Once this is done, no further action is needed! The simulation objects will automatically be moved to the GPU and back once the simulation is finished. When the simulation is run a message will be shown with either the GPU device being used or the number of CPU threads if running on the CPU.

Of course, it is still possible to move objects to the GPU manually, and control precision using the f32 and f64 functions:

```julia
x = rand(100)
x |> f32 |> gpu # Float32 CuArray
```

To change the precision level used for the entire simulation, the `sim_params["precision"]` parameter can be set to either `f32` or `f64` (Note that for most GPUs, Float32 operations are considerably faster compared with Float64). In addition, the `sim_params["gpu"]` option can be set to true or false to enable / disable the gpu functionality (if set to true, the backend package will still need to be loaded beforehand):

```julia
using KomaMRI
using CUDA
sys = Scanner
obj = brain_phantom2D()
seq = PulseDesigner.EPI_example()

#Simulate on the GPU using 32-bit floating point values
sim_params = Dict{String,Any}(
  "Nblocks" => 20,
  "gpu" => true,
  "precision" => "f32"
  "sim_method" => Bloch(),
)
simulate(obj, seq, sys; sim_params)
```


## How Objects are moved to the GPU

Koma's `gpu` function implementation calls a separate `gpu` function with a backend parameter of type `<:KernelAbstractions.GPU` for the backend it is using. This function then calls the `fmap` function from package `Functors.jl` to recursively call `adapt` from package `Adapt.jl` on each field of the object being transferred. This is similar to how many other Julia packages, such as `Flux.jl`, transfer data to the GPU. However, an important difference is that KomaMRI adapts directly to the `KernelAbstractions.Backend` type in order to use the `adapt_storage` functions defined in each backend package, rather than defining custom adapters, resulting in an implementation with fewer lines of code.

## Inside the Simulation

KomaMRI has three different simulation methods, all of which can run on the GPU: 

* `BlochSimple`: [BlochSimple.jl](https://github.com/JuliaHealth/KomaMRI.jl/blob/master/KomaMRICore/src/simulation/SimMethods/BlochSimple/BlochSimple.jl)
* `BlochDict`: [BlochDict.jl](https://github.com/JuliaHealth/KomaMRI.jl/blob/master/KomaMRICore/src/simulation/SimMethods/BlochDict/BlochDict.jl)
* `Bloch`: [BlochCPU.jl](https://github.com/JuliaHealth/KomaMRI.jl/blob/master/KomaMRICore/src/simulation/SimMethods/Bloch/BlochCPU.jl) / [BlochGPU.jl](https://github.com/JuliaHealth/KomaMRI.jl/blob/master/KomaMRICore/src/simulation/SimMethods/Bloch/BlochGPU.jl)

`BlochSimple` is the simplest method and prioritizes readability. 

`BlochDict` can be understood as an extension to `BlochSimple` that outputs a more detailed signal.

`Bloch` is equivalent to `BlochSimple` in the operations it performs, but is much faster since it has been optimized both for the CPU and GPU. The CPU implementation prioritizes conserving memory, and makes extensive use of pre-allocation for the simulation arrays. Unlike the GPU implementation, it does not allocate a matrix of size `Number of Spins x Number of Time Points` in each block, instead using a for loop to step through time.

In contrast, the GPU implementation divides work among as many threads as possible at the beginning of the `run_spin_precession!` and `run_spin_excitation!` functions. For the CPU implementation, this would not be beneficial since there are far less CPU threads available compared with the GPU. Preallocation is also used via the same `prealloc` function used in `BlochCPU.jl`, where a struct of arrays is allocated at the beginning of the simulation that can be re-used in each simulation block. In addition, a `precalc` function is called before moving the simulation objects to the GPU to do certain calculations that are faster on the CPU beforehand.

Compared with `BlochSimple`, which only uses array broadcasting for parallelization, `Bloch` also uses kernel-based methods in its `run_spin_excitation!` function for operations which need to be done sequentially. The [kernel implementation](https://github.com/JuliaHealth/KomaMRI.jl/blob/master/KomaMRICore/src/simulation/SimMethods/Bloch/KernelFunctions.jl) uses shared memory to store the necessary arrays for applying the spin excitation for fast memory access, and separates the complex arrays into real and imaginary components to avoid bank conflicts.

The performance differences between Bloch and BlochSimple can be seen on the KomaMRI [benchmarks page](https://juliahealth.org/KomaMRI.jl/benchmarks/). The first data point is from when `Bloch` was what is now `BlochSimple`, before a more optimized implementation was created. The following three pull requests are primarily responsible for the performance differences between `Bloch` and `BlochSimple`:

* [(443) Optimize run_spin_precession! and run_spin_excitation! for CPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/443)
* [(459) Optimize run_spin_precession! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/459)
* [(462) Optimize run_spin_excitation! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/462)
