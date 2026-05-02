function progressbar_callback(Nblocks)
    progress_bar = Progress(Nblocks; desc="Running simulation...")
    function progressbar_callback_fun(progress_info, simulation_blocks_info, device_data, sim_params)
        next!(
            progress_bar;
            showvalues=[
                (:simulated_blocks, progress_info.block), 
                (:rf_blocks, progress_info.rfs), 
                (:acq_samples, progress_info.samples - 1)
            ],
        )
    end
    return progressbar_callback_fun
end