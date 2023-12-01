"""
Update the encoding parameters of the RawAcquisitionData variable
In the sequence must be present the "Nx" parameter
"""
function update_encoding_limits!(raw::RawAcquisitionData, seq::Sequence)
    Nx = get(seq.DEF, "Nx", 1)
    encoding_limits = Dict(
        "enc_lim_kspace_encoding_step_1" => Limit(0, Nx-1, ceil(Int, Nx / 2)),   #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_2" => Limit(0, 0, 0),     #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, 0, 0),     #min, max, center, e.g. imaging slice number
        "enc_lim_contrast"               => Limit(0, 0, 0),     #min, max, center, e.g. echo number in multi-echo
        "enc_lim_phase"                  => Limit(0, 0, 0),     #min, max, center, e.g. cardiac phase number
        "enc_lim_repetition"             => Limit(0, 0, 0),     #min, max, center, e.g. dynamic number for dynamic scanning
        "enc_lim_set"                    => Limit(0, 0, 0),     #min, max, center, e.g. flow encoding set
        "enc_lim_segment"                => Limit(0, 0, 0),     #min, max, center, #segment: e.g. segment number for segmented acquisition
    )
    merge!(raw.params, encoding_limits)
end
