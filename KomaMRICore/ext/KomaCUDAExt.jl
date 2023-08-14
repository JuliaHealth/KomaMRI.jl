module KomaCUDAExt

    using KomaMRICore
    using KomaMRICore: KomaCPUAdaptor, KomaCUDAAdaptor, fmap
    using CUDA
    using Adapt
    import KomaMRICore: _cuda       #To extend it
    import Adapt: adapt_storage

    const USE_CUDA = Ref{Union{Nothing, Bool}}(nothing)

    function check_use_cuda()
        if !isnothing(USE_CUDA[])
            return
        end
        USE_CUDA[] = CUDA.functional()
        if !USE_CUDA[]
            @info """
            The CUDA function is being called but CUDA.jl is not functional.
            Defaulting back to the CPU. (No action is required if you want to run on the CPU).
            """ maxlog=1
        end
        return
    end

    function __init__()
        KomaMRICore.CUDA_LOADED[] = true
    end

    # GPU adaptor
    adapt_storage(to::KomaCUDAAdaptor, x) = CUDA.cu(x)

    function _cuda(x)
        check_use_cuda()
        USE_CUDA[] || return x
        fmap(x -> Adapt.adapt(KomaCUDAAdaptor(), x), x; exclude=KomaMRICore._isleaf)
    end

    """
    print_gpus()

    Simple function to print the CUDA devices available in the host.
    """
    print_gpus() = begin
        check_use_cuda()
        if USE_CUDA[]
            println( "$(length(devices())) CUDA capable device(s)." )
            for (i,d) = enumerate(devices())
                u = i == 1 ? "*" : " "
                println( "  ($(i-1)$u) $(name(d))")
            end
        else
            println("0 CUDA capable devices(s).")
        end
    end

end