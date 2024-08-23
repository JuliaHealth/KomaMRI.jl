# GPU Parallelization

KomaMRI uses a vendor agnostic approach to GPU parallelization in order to support multiple GPU backends. Currently, the following backends are supported:

* CUDA.jl (Nvidia)
* Metal.jl (Apple)
* AMDGPU.jl (AMD)
* oneAPI.jl (Intel)

## Choosing a GPU Backend

To determine which backend to use, KomaMRI uses [package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) (introduced in Julia 1.9) to avoid having the packagaes for each GPU backend as explicit dependencies. This means that the user is responsible for loading the backend package (e.g. `using CUDA`) at the beginning of their code, or prior to calling KomaUI(), otherwise, Koma will default back to the CPU. Once this is done, no further action is recquired! The simulation objects will automatically be moved to the GPU and back once the simulation is finished. When the simulation is run a message will be shown with either the GPU device being used or the number of CPU threads if running on the CPU.

## How Objects are moved to the GPU

KomaMRI has a general purpose function, `gpu`, to move data from the CPU to the GPU. The `gpu` function implementation calls a separate `gpu` function with a backend parameter of type `<:KernelAbstractions.GPU` for the backend it is using. This function then calls the `fmap` function from package `Functors.jl` to recursively call `adapt` from package `Adapt.jl` on each field of the object being transferred. This is similar to how many other Julia packages, such as `Flux.jl`, transfer data to the GPU. However, an important difference is that KomaMRI adapts directly to the `KernelAbstractions.Backend` type in order to use the `adapt_storage` functions defined in each backend package, rather than defining custom adapters, resulting in an implementation with fewer lines of code.

## Inside the Simulation

KomaMRI has three different simulation methods, all of which can run on the GPU: 

* `Bloch`
* `BlochSimple`
* `BlochDict`

Of the three methods, `Bloch` is the most optimized, and has separate implementations specialized for the CPU and GPU. `BlochSimple` is equivalent to `Bloch` in the operations it performs, but less optimized and easier to understand. `BlochDict` can be understood as an extension of `BlochSimple` that outputs a more complete signal. 

`BlochSimple` and `Bloch` take slightly different approaches to GPU parallelization. `BlochSimple` exclusively uses array broadcasting, with parallelization on the arrays being done implicitly by the GPU compiler. In constrast, `Bloch` uses explicit GPU kernels where advantageous, using package `KernelAbstractions.jl`. Readers curious about the performance improvements between `Bloch` and `BlochSimple` may want to look at the following pull reqeusts:

* [(459) Optimize run_spin_precession! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/459)
* [(462) Optimize run_spin_excitation! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/462)
