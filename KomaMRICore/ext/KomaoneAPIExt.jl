module KomaoneAPIExt

using oneAPI
import KomaMRICore
import Adapt

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name
@inline KomaMRICore._cis(x) = cos(x) + im * sin(x)

function KomaMRICore._print_devices(::oneAPIBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => oneAPI.properties(d).name for
        (i, d) in enumerate(oneAPI.devices())
    ]
    @info "$(length(oneAPI.devices())) oneAPI capable device(s)." devices...
end

#Temporary workaround since oneAPI.jl (similar to Metal) does not support some array operations
#Once run_spin_excitation! and run_spin_precession! are kernel-based, this code can be removed
Base.cumsum(x::oneVector{T}) where T = convert(oneVector{T}, cumsum(KomaMRICore.cpu(x)))
Base.findall(x::oneVector{Bool}) = convert(oneVector, findall(KomaMRICore.cpu(x)))

using KernelAbstractions: @index, @kernel, @Const, synchronize

"""Naive cumsum implementation for matrix, parallelizes along the first dimension"""
Base.cumsum(A::oneArray{T}; dims) where T = begin
    dims == 2 || @error "oneAPI cumsum implementation only supports keyword argument dims=2"
    backend = oneAPIBackend()
    B = similar(A)
    cumsum_kernel = naive_cumsum!(backend)
    cumsum_kernel(B, A, ndrange=size(A,1))
    synchronize(backend)
    return B
end

Base.cumsum!(B::oneArray{T}, A::oneArray{T}; dims) where T = begin
    dims == 2 || @error "oneAPI cumsum implementation only supports keyword argument dims=2"
    backend = oneAPIBackend()
    cumsum_kernel = naive_cumsum!(backend)
    cumsum_kernel(B, A, ndrange=size(A,1))
    synchronize(backend)
end

## COV_EXCL_START
@kernel function naive_cumsum!(B, @Const(A))
    i = @index(Global)
    T = eltype(A)

    cur_val = zero(T)
    for k ∈ 1:size(A, 2)
        @inbounds cur_val += A[i, k]
        @inbounds B[i, k] = cur_val
    end
end
## COV_EXCL_STOP

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], oneAPIBackend())
    @warn "oneAPI does not support all array operations used by KomaMRI. GPU performance may be slower than expected"
end

end