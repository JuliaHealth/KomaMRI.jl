const LOADED_BACKENDS = Ref{Vector{KA.GPU}}([])

name(::KA.CPU) = "CPU"
isfunctional(::KA.CPU) = true
isfunctional(x) = false
set_device!(backend, val) = @error "set_device! called with invalid parameter types: '$(typeof(backend))', '$typeof(val)'" 
gpu_name(backend) = @error "gpu_name called with invalid backend type $(typeof(backend))"
function reclaim_gpu(backend) end

function get_backend()
    if isempty(LOADED_BACKENDS[])
        @info """ 
        The GPU functionality is being called but a GPU backend must be loaded
        to access it. Add 'using CUDA / Metal / AMDGPU / oneAPI' to your code.
        Defaulting back to the CPU. (No action is required if you want to run on the CPU).
        """ maxlog=1
        return KA.CPU()
    end

    functional_backends = []
    for backend in LOADED_BACKENDS[]
        if isfunctional(backend)
            push!(functional_backends, backend)
        else
            @warn "Loaded backend $(name(backend)) is not functional"
        end
    end
    
    if length(functional_backends) == 1
        @info """Using  backend: '$name(backend)'"""  maxlog = 1
        return functional_backends[1]
    elseif length(functional_backends) == 0
        @info """ Defaulting back to the CPU. (No action is required if you want to run on the CPU). """ maxlog = 1
        return KA.CPU()
    else
        # Will probably never get here
        @info """
        Multiple functional backends have been loaded and KomaMRI does not know which one
        to use. Ensure that your code contains only one 'using' statement for the GPU backend
        you wish to use. Defaulting back to the CPU. (No action is required if you want to run 
        on the CPU).
        """ maxlog = 1
        return KA.CPU()
    end
end

