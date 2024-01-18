import Functors: @functor, functor, fmap, isleaf
import Adapt: adapt, adapt_storage
#Checks if CUDA is available for the session
const use_cuda = Ref{Union{Nothing,Bool}}(nothing)

"""
    print_gpus()

Simple function to print the CUDA devices available in the host.
"""
function print_gpus()
    check_use_cuda()
    if use_cuda[]
	    cuda_devices = [Symbol("($(i-1)$(i == 1 ? "*" : " "))") => name(d) for (i,d) = enumerate(devices())]
		@info "$(length(devices())) CUDA capable device(s)." cuda_devices...
    else
        @info "0 CUDA capable devices(s)."
    end
	return
end

"""
Checks if the PC has a functional CUDA installation. Inspired by Flux's `check_use_cuda` funciton.
"""
function check_use_cuda()
	if use_cuda[] === nothing
		use_cuda[] = CUDA.functional()
		if !(use_cuda[])
		@info """The GPU function is being called but the GPU is not accessible.
					Defaulting back to the CPU. (No action is required if you want to run on the CPU).""" maxlog=1
		end
	end
end

#Aux. funcitons to check if the variable we want to convert to CuArray is numeric
_isbitsarray(::AbstractArray{<:Real}) = true
_isbitsarray(::AbstractArray{T}) where T = isbitstype(T)
_isbitsarray(x) = false
_isleaf(x) = _isbitsarray(x) || isleaf(x)

# GPU adaptor
struct KomaCUDAAdaptor end
adapt_storage(to::KomaCUDAAdaptor, x) = CUDA.cu(x)

"""
	gpu(x)

Tries to move `x` to the current GPU device. Inspired by Flux's `gpu` function.

This works for functions, and any struct marked with `@functor`.

Use [`cpu`](@ref) to copy back to ordinary `Array`s.

See also [`f32`](@ref) and [`f64`](@ref) to change element type only.

# Examples
```julia
x = x |> gpu
```
"""
function gpu(x)
	check_use_cuda()
	use_cuda[] ? fmap(x -> adapt(KomaCUDAAdaptor(), x), x; exclude = _isleaf) : x
end

#CPU adaptor
struct KomaCPUAdaptor end
adapt_storage(to::KomaCPUAdaptor, x::AbstractArray) = adapt(Array, x)
adapt_storage(to::KomaCPUAdaptor, x::AbstractRange) = x

# To CPU
"""
	cpu(x)

Tries to move object to CPU. Inspired by Flux's `cpu` function.

This works for functions, and any struct marked with `@functor`.

See also [`gpu`](@ref).

# Examples
```julia
x = x |> cpu
```
"""
cpu(x) = fmap(x -> adapt(KomaCPUAdaptor(), x), x)

#Precision
paramtype(T::Type{<:Real}, m) = fmap(x -> adapt(T, x), m)

adapt_storage(T::Type{<:Real}, xs::AbstractArray{<:Real}) = convert.(T, xs) #Type piracy
adapt_storage(T::Type{<:Real}, xs::AbstractArray{<:Complex}) = convert.(Complex{T}, xs) #Type piracy
adapt_storage(T::Type{<:Real}, xs::AbstractArray{<:Bool}) = xs #Type piracy

"""
    f32(m)

Converts the `eltype` of model's parameters to `Float32`.
Recurses into structs marked with `@functor`.

See also [`f64`](@ref).
"""
f32(m) = paramtype(Float32, m)

"""
    f64(m)

Converts the `eltype` of model's parameters to `Float64`.
Recurses into structs marked with `@functor`.

See also [`f32`](@ref).
"""
f64(m) = paramtype(Float64, m)

#The functor macro makes it easier to call a function in all the parameters
@functor Phantom
@functor Spinor
@functor DiscreteSequence

#Exporting functions
export gpu, cpu, f32, f64
