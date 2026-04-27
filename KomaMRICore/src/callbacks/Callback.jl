struct Callback{F}
    enabled::Bool
    every::Int
    fun::F
end

Callback(every::Int, fun::F) where {F} = Callback(every > 0, every, fun)

function (cb::Callback)(progress_info, simulation_blocks_info, device_data, sim_params)
    if cb.enabled && (progress_info.block - 1) % cb.every == 0
        cb.fun(progress_info, simulation_blocks_info, device_data, sim_params)
    end
end

# Included Callbacks
include("progress_callback.jl")