function ui_progressbar_callback(w)
    function ui_progressbar_callback_fun(progress_info, simulation_blocks_info, device_data, sim_params)
        update_blink_window_progress!(w, progress_info.block, progress_info.Nblocks)
    end
    return Callback(1, ui_progressbar_callback_fun)
end
