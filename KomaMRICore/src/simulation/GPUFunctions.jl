import Functors: @functor, functor, fmap, isleaf
import Adapt: adapt, adapt_storage
using Preferences

# CPU adaptor
struct KomaCPUAdaptor end

# Define rules for handling structured arrays
adapt_storage(to::KomaCPUAdaptor, x::AbstractArray) = adapt(Array, x)
adapt_storage(to::KomaCPUAdaptor, x::AbstractRange) = x

# To CPU
"""
	cpu(x)

Tries to move object to CPU. Inspired by Koma's `cpu` function.

This works for functions, and any struct marked with [`@functor`](@ref).

# Examples
```julia
x = x |> cpu
```
"""
cpu(x) = fmap(x -> adapt(KomaCPUAdaptor(), x), x)

# Aux. funcitons to check if the variable we want to convert to CuArray is numeric
_isbitsarray(::AbstractArray{<:Real}) = true
_isbitsarray(::AbstractArray{T}) where T = isbitstype(T)
_isbitsarray(x) = false
_isleaf(x) = _isbitsarray(x) || isleaf(x)

# Available backends for GPU acceleration
const GPU_BACKENDS = ("CUDA", "AMD", "Metal")
const GPU_BACKEND = @load_preference("gpu_backend", "CUDA")

# Choose GPU backend
function gpu_backend!(backend::String)
    if backend == GPU_BACKEND
        @info """
        GPU backend is already set to: $backend.
        No need to do anything else.
        """
        return
    end
    backend in GPU_BACKENDS || throw(ArgumentError("""
    Unsupported GPU backend: $backend.
    Supported backends are: $GPU_BACKENDS.
    """))
    @set_preferences!("gpu_backend" => backend)
    @info """
    New GPU backend set: $backend.
    Restart your Julia session for this change to take effect!
    """
end

"""
	gpu(x)

Tries to move `x` to the current GPU device. Inspired by Koma's `gpu` function.

This works for functions, and any struct marked with [`@functor`](@ref).

# Examples
```julia
x = x |> gpu
```
"""
function gpu(x)
    @static if GPU_BACKEND == "CUDA"
        gpu(KomaCUDAAdaptor(), x)
    elseif GPU_BACKEND == "AMD"
        gpu(KomaAMDAdaptor(), x)
    elseif GPU_BACKEND == "Metal"
        gpu(KomaMetalAdaptor(), x)
    else
        error("""
        Unsupported GPU backend: $GPU_BACKEND.
        Supported backends are: $GPU_BACKENDS.
        """)
    end
end

#Precision adaptor
struct KomaEltypeAdaptor{T} end

# Define rules for handling structured arrays
adapt_storage(::KomaEltypeAdaptor{T}, x::AbstractArray{<:AbstractFloat}) where {T<:AbstractFloat} = convert(AbstractArray{T}, x)
adapt_storage(::KomaEltypeAdaptor{T}, x::AbstractArray{<:Complex{<:AbstractFloat}}) where {T<:AbstractFloat} = convert(AbstractArray{Complex{T}}, x)
# adapt_storage(::KomaEltypeAdaptor{T}, x::AbstractArray{<:Bool}) where {T<:AbstractFloat} = x

_paramtype(::Type{T}, x) where T = fmap(adapt(KomaEltypeAdaptor{T}()), x)

# Fastpath for arrays
_paramtype(::Type{T}, x::AbstractArray{<:AbstractFloat}) where {T<:AbstractFloat} = convert(AbstractArray{T}, x)
_paramtype(::Type{T}, x::AbstractArray{<:Complex{<:AbstractFloat}}) where {T<:AbstractFloat} = convert(AbstractArray{Complex{T}}, x)

"""
    f32(x)
Converts the `eltype` of model's parameters to `Float32`
Recurses into structs marked with [`@functor`](@ref).
"""
f32(x) = _paramtype(Float32, x)

"""
    f64(x)
Converts the `eltype` of model's parameters to `Float64` (which is Koma's default)..
Recurses into structs marked with [`@functor`](@ref).
"""
f64(x) = _paramtype(Float64, x)

######## GPU MULTIVENDOR SUPPORT ########
######## CUDA extension
struct KomaCUDAAdaptor end
const CUDA_LOADED = Ref{Bool}(false)
function gpu(::KomaCUDAAdaptor, x)
    if CUDA_LOADED[]
        return _cuda(x)
    else
        @info """
        The CUDA functionality is being called but
        `CUDA.jl` must be loaded to access it.
        Add `using CUDA` or `import CUDA` to your code.
        """ maxlog=1
        return x
    end
end
function _cuda end

######## AMDGPU extension
struct KomaAMDAdaptor end
const AMDGPU_LOADED = Ref{Bool}(false)
function gpu(::KomaAMDAdaptor, x)
    if AMDGPU_LOADED[]
        return _amd(x)
    else
        @info """
        The AMDGPU functionality is being called but
        `AMDGPU.jl` must be loaded to access it.
        Add `using AMDGPU` or `import AMDGPU` to your code.
        """ maxlog=1
        return x
    end
end
function _amd end

######## Metal extension. Apple M1/M2  
struct KomaMetalAdaptor end
const METAL_LOADED = Ref{Bool}(false)
function gpu(::KomaMetalAdaptor, x)
    if METAL_LOADED[]
        return _metal(x)
    else
        @info """
        The Metal functionality is being called but
        `Metal.jl` must be loaded to access it.
        """ maxlog=1
        return x
    end
end
function _metal end

"""
print_gpus()

Simple function to print the CUDA devices available in the host.
"""
function print_gpus()
    println("No GPU_BACKEND has been loaded.")
end

#The functor macro makes it easier to call a function in all the parameters
@functor Phantom
@functor Spinor
@functor DiscreteSequence

#Exporting functions_paramtype
export gpu, cpu, f32, f64
