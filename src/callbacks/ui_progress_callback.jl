function ui_progressbar_callback(w)
    return Callback(1, (progress_info, _, _, _) -> begin
        update_blink_window_progress!(w, progress_info.block, progress_info.Nblocks)
    end)
end
