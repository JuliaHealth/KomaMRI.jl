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

Two other simulation parameters, `gpu_groupsize_precession` and `gpu_groupsize_excitation` are exposed to allow adjusting the number of threads in each threadgroup within the `run_spin_precession!` and `run_spin_excitation!` gpu kernels. By default, they are both 256, however, on some devices other values may result in faster performance. The gpu groupsize must be a multiple of 32 between 32 and no higher than 1024, otherwise the simulation will throw an error.

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

`BlochSimple` uses array-based methods which are simpler to understand compared with the more optimized `Bloch` implementation.

`BlochDict` can be understood as an extension to `BlochSimple` that outputs a more detailed signal.

`Bloch` is equivalent to `BlochSimple` in the operations it performs, but has separate implementations optimized for both the CPU and GPU. The CPU implementation uses array broadcasting for computation and preallocates all simulation arrays to conserve memory. The simulation arrays are 1-dimensional with length equal to the number of spins in the phantom and are updated at each time step. The GPU implementation also uses preallocation and a similar loop-based computation strategy, but does so using kernels for spin precession and excitation implemented using the `KernelAbstractions.jl` package. A key advantage of using kernel-based methods is that intermediate values compuated based on phantom and sequence properties can be stored in registers without having to write back to GPU global memory, which has much higher memory latency compared with the CPU. Other optimizations within the kernels include:

* Reducing the output signal value at each time step within the kernel so that the first thread for each thread group writes the sum of the signal values for each thread in the threadgroup to GPU global memory. This reduces the number of GPU global memory reads + writes needed for the output signal from `Number of Spins x Number of Time Points` to `Number of Spins x Number of Time Points / Number of Threads in Threadgroup`, improving scalability for large phantom objects.

* Using julia's `Val` type to specialize at compile-time on properties unique to the simulation inputs. For example, whether the phantom exibits spin motion is passed as either `Val(true)` or `Val(false)` to the precession kernel so that different kernels will be compiled for phantoms with or without motion. For the kernel compiled for phantoms without motion, there will be no runtime check of the motion type of the phantom, and everything inside the `if MOTION` statements in the kernel will be compiled out, saving register space and enabling further compiler optimizations. This strategy enables adding support in the future for less common use cases without negatively impacting performance for simulations not using these features.

* Since GPU registers are limited and can hurt GPU occupancy if a kernel uses a high number, their use is minimized by working with real and imaginary components directly rather than abstracting complex number math, using unsigned int32 literal values instead of Julia's default `Int64`, and inlining all functions called from within the kernels.

The performance differences between Bloch and BlochSimple can be seen on the KomaMRI [benchmarks page](https://juliahealth.org/KomaMRI.jl/benchmarks/). The first data point is from when `Bloch` was what is now `BlochSimple`, before a more optimized implementation was created. The following pull requests are primarily responsible for the performance differences between `Bloch` and `BlochSimple`:

* [(443) Optimize run_spin_precession! and run_spin_excitation! for CPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/443)
* [(459) Optimize run_spin_precession! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/459)
* [(462) Optimize run_spin_excitation! for GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/462)
* [(537) Faster Bloch GPU](https://github.com/JuliaHealth/KomaMRI.jl/pull/537)
