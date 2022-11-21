using CUDA

const dim = 1_000_000
device!(1)
x = CUDA.ones(Float32, dim)
y = CUDA.ones(Float32, dim)
z = CUDA.zeros(Float32, dim)

function saxp_gpu_kernel!(z, x, y)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= length(z)
        @inbounds z[i] = x[i] + y[i]
    end
    return nothing
end

nthreads = CUDA.attribute(
    device(),
    CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK
)
nblocks = cld(dim, nthreads) #Ceiling division

@time CUDA.@sync @cuda(
    threads = nthreads,
    blocks = nblocks,
    saxp_gpu_kernel!(z, x, y)
)

@time CUDA.@sync z .= x .+ y;