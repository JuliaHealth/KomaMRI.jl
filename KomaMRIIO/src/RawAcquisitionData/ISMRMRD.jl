"""
    raw = signal_to_raw_data(signal, seq; phantom_name, sys, sim_params)

Transforms the raw signal into a RawAcquisitionData struct (nearly equivalent to the ISMRMRD
format) used for reconstruction with MRIReco.

# Arguments
- `signal`: (`::Matrix{Complex}`) raw signal matrix
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `phantom_name`: (`::String`, `="Phantom"`) phantom name
- `sys`: (`::Scanner`, `=Scanner()`) Scanner struct
- `sim_params`: (`::Dict{String, Any}`, `=Dict{String,Any}()`) simulation parameter dictionary

# Returns
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_se.seq")

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> sim_params = KomaMRICore.default_sim_params(); sim_params["return_type"] = "mat"

julia> signal = simulate(obj, seq, sys; sim_params)

julia> raw = signal_to_raw_data(signal, seq)

julia> plot_signal(raw)
```
"""
function signal_to_raw_data(
    signal, seq;
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}()
)
    version = string(VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["version"]))
    #Number of samples and FOV
    _, ktraj = get_kspace(seq) #kspace information
    mink = minimum(ktraj, dims=1)
    maxk = maximum(ktraj, dims=1)
    Wk = maxk .- mink
    Δx = 1 ./ Wk[1:2] #[m] Only x-y
    Nx = get(seq.DEF, "Nx", 1)
    Ny = get(seq.DEF, "Ny", 1)
    Nz = get(seq.DEF, "Nz", 1)
    if haskey(seq.DEF, "FOV")
        FOVx, FOVy, _ = seq.DEF["FOV"] #[m]
        if FOVx > 1 FOVx *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        if FOVy > 1 FOVy *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        Nx = round(Int64, FOVx / Δx[1])
        Ny = round(Int64, FOVy / Δx[2])
    else
        FOVx = Nx * Δx[1]
        FOVy = Ny * Δx[2]
    end
    #It needs to be transposed for the raw data
    ktraj = maximum(2*abs.(ktraj[:])) == 0 ? transpose(ktraj) : transpose(ktraj)./ maximum(2*abs.(ktraj[:]))

    #First we define the ISMRMRD data XML header
    params = Dict(
        #Simulator
        "systemVendor"                   => "KomaMRI.jl", #String
        "systemModel"                    => "v"*version, #String
        "institutionName"                => "Pontificia Universidad Catolica de Chile", #String
        #Phantom
        "patientName"                    => phantom_name,
        #Sys
        "systemFieldStrength_T"          => sys.B0, #Float
        "H1resonanceFrequency_Hz"        => floor(Int64, γ * sys.B0), #Long (Int)
        #Seq
        "trajectory"                     => "other", #Must be: cartesian, epi, radial, goldenangle, spiral, and other
        "protocolName"                   => haskey(seq.DEF,"Name") ? seq.DEF["Name"] : "NoName", #String
        # "trajectoryDescription"          => Dict{String, Any}("comment"=>""), #You can put wathever you want here: comment, bandwidth, MaxGradient_G_per_cm, MaxSlewRate_G_per_cm_per_s, interleaves, etc
        "userParameters"                 => sim_params, #Dict with parameters
        #encoding>
        "encodedFOV"                     => [FOVx, FOVy, 1.0],    #encodedSpace>fieldOfView_mm
        "reconFOV"                       => [FOVx, FOVy, 1.0],    #reconSpace>fieldOfView_mm
        "encodedSize"                    => [Nx, Ny, 1],            #encodedSpace>matrixSize
        "reconSize"                      => [Nx+Nx%2, Ny+Ny%2, 1],  #reconSpace>matrixSize
        #encodingLimits>
        "enc_lim_kspace_encoding_step_1" => Limit(0, Nx-1, ceil(Int, Nx / 2)),   #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_2" => Limit(0, 0, 0),     #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, 0, 0),     #min, max, center, e.g. imaging slice number
        "enc_lim_contrast"               => Limit(0, 0, 0),     #min, max, center, e.g. echo number in multi-echo
        "enc_lim_phase"                  => Limit(0, 0, 0),     #min, max, center, e.g. cardiac phase number
        "enc_lim_repetition"             => Limit(0, 0, 0),     #min, max, center, e.g. dynamic number for dynamic scanning
        "enc_lim_set"                    => Limit(0, 0, 0),     #min, max, center, e.g. flow encoding set
        "enc_lim_segment"                => Limit(0, 0, 0),     #min, max, center, #segment: e.g. segment number for segmented acquisition
        #sequenceParameters
        # "TR"                             => 0,
        # "TE"                             => 0,
        # "TI"                             => 0,
        # "flipAngle_deg"                  => 0,
        # "echo_spacing"                   => 0,
    )

    #Then, we define the Profiles
    profiles = Profile[]
    t_acq = KomaMRICore.get_adc_sampling_times(seq)
    Nadcs = sum(is_ADC_on.(seq))
    NadcsPerImage = floor(Int, Nadcs / Nz)
    scan_counter = 0
    nz = 0
    current = 1
    for s = seq #Iterate over sequence blocks
        if is_ADC_on(s)
            Nsamples = s.ADC.N[1]
            Δt_us = Float32( s.ADC.T[1] / (Nsamples - 1) * 1e6 )
            t0_us = floor(Int32, t_acq[current] * 1e6 )
            #Header of profile data, head::AcquisitionHeader
            head = AcquisitionHeader(
                0, #version: First unsigned int indicates the version
                0, #flags: bit field with flags, noise measure, calibration, coil sens, measure, etc
                0, #measurement_uid: Unique ID for the measurement
                scan_counter, #scan_counter: Current acquisition number in the measurement
                t0_us, #acquisition_time_stamp: Acquisition clock, I am "miss"-using this variable to store t0 in us
                (0, 0, 0), #physiology_time_stamp: Physiology time stamps, e.g. ecg, breating, etc.
                Nsamples, #number_of_samples
                1, #available_channels: Available coils
                1, #active_channels: Active coils on current acquisiton
                Tuple(0 for i=1:16), #channel_mask: Active coils on current acquisiton
                0, #discard_pre: Samples to be discarded at the beginning of acquisition
                0, #discard_post: Samples to be discarded at the end of acquisition
                0, #center_sample: Sample at the center of k-space
                0, #encoding_space_ref: Reference to an encoding space, typically only one per acquisition
                2, #trajectory_dimensions: Indicates the dimensionality of the trajectory vector (0 means no trajectory)
                Δt_us, #sample_time_us: Time between samples in micro seconds, sampling BW
                (0.0f0, 0.0f0, 0.0f0), #position: Three-dimensional spatial offsets from isocenter
                (1.0f0, 0.0f0, 0.0f0), #read_dir: Directional cosines of the readout/frequency encoding
                (0.0f0, 1.0f0, 0.0f0), #phase_dir: Directional cosines of the phase
                (0.0f0, 0.0f0, 1.0f0), #slice_dir: Directional cosines of the slice direction
                (0.0f0, 0.0f0, 0.0f0), #patient_table_position: Patient table off-center
                EncodingCounters( #idx: Encoding loop counters
                    scan_counter, #kspace_encode_step_1: e.g. phase encoding line number
                    0, #kspace_encode_step_2: e.g. partition encoding number
                    0, #average: e.g. signal average number
                    nz, #slice: e.g. imaging slice number
                    0, #contrast: e.g. echo number in multi-echo
                    0, #phase: e.g. cardiac phase number
                    0, #repetition: e.g. dynamic number for dynamic scanning
                    0, #set: e.g. flow encoding set
                    0, #segment: e.g. segment number for segmented acquisition
                    (0, 0, 0, 0, 0, 0, 0, 0) #user: Free user parameters
                ),
                (0, 0, 0, 0, 0, 0, 0, 0), #user_int: Free user parameters
                (0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0) #user_float: Free user parameters
            )
            #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples
            traj = ktraj[1:2, current:current+Nsamples-1]
            #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
            dat =  signal[current:current+Nsamples-1, :]
            #Saving profile
            push!(profiles, Profile(head, traj, dat))
            #Update counters
            scan_counter += 1
            current += Nsamples
            if scan_counter % NadcsPerImage == 0 #For now only Nz is considered
                nz += 1 #another image
                scan_counter = 0 #reset counter
            end
        end
    end

    return RawAcquisitionData(params, profiles)
end
# enum ISMRMRD_AcquisitionFlags {
#     ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1               =  1,
#     ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1                =  2,
#     ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2               =  3,
#     ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2                =  4,
#     ISMRMRD_ACQ_FIRST_IN_AVERAGE                    =  5,
#     ISMRMRD_ACQ_LAST_IN_AVERAGE                     =  6,
#     ISMRMRD_ACQ_FIRST_IN_SLICE                      =  7,
#     ISMRMRD_ACQ_LAST_IN_SLICE                       =  8,
#     ISMRMRD_ACQ_FIRST_IN_CONTRAST                   =  9,
#     ISMRMRD_ACQ_LAST_IN_CONTRAST                    = 10,
#     ISMRMRD_ACQ_FIRST_IN_PHASE                      = 11,
#     ISMRMRD_ACQ_LAST_IN_PHASE                       = 12,
#     ISMRMRD_ACQ_FIRST_IN_REPETITION                 = 13,
#     ISMRMRD_ACQ_LAST_IN_REPETITION                  = 14,
#     ISMRMRD_ACQ_FIRST_IN_SET                        = 15,
#     ISMRMRD_ACQ_LAST_IN_SET                         = 16,
#     ISMRMRD_ACQ_FIRST_IN_SEGMENT                    = 17,
#     ISMRMRD_ACQ_LAST_IN_SEGMENT                     = 18,
#     ISMRMRD_ACQ_IS_NOISE_MEASUREMENT                = 19,
#     ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION             = 20,
#     ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 21,
#     ISMRMRD_ACQ_IS_REVERSE                          = 22,
#     ISMRMRD_ACQ_IS_NAVIGATION_DATA                  = 23,
#     ISMRMRD_ACQ_IS_PHASECORR_DATA                   = 24,
#     ISMRMRD_ACQ_LAST_IN_MEASUREMENT                 = 25,
#     ISMRMRD_ACQ_IS_HPFEEDBACK_DATA                  = 26,
#     ISMRMRD_ACQ_IS_DUMMYSCAN_DATA                   = 27,
#     ISMRMRD_ACQ_IS_RTFEEDBACK_DATA                  = 28,
#     ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 29,
#     ISMRMRD_ACQ_COMPRESSION1                        = 53,
#     ISMRMRD_ACQ_COMPRESSION2                        = 54,
#     ISMRMRD_ACQ_COMPRESSION3                        = 55,
#     ISMRMRD_ACQ_COMPRESSION4                        = 56,
#     ISMRMRD_ACQ_USER1                               = 57,
#     ISMRMRD_ACQ_USER2                               = 58,
#     ISMRMRD_ACQ_USER3                               = 59,
#     ISMRMRD_ACQ_USER4                               = 60,
#     ISMRMRD_ACQ_USER5                               = 61,
#     ISMRMRD_ACQ_USER6                               = 62,
#     ISMRMRD_ACQ_USER7                               = 63,
#     ISMRMRD_ACQ_USER8                               = 64
# };

Base.show(io::IO, raw::RawAcquisitionData) = begin
    Nt, Nc = size(raw.profiles[1].data)
    compact = get(io, :compact, false)
    seq_name = get(raw.params, "protocolName", "None")
	if !compact
        print(io, "RawAcquisitionData[SeqName: $seq_name | $(length(raw.profiles)) Profile(s) of $Nt×$Nc]")
    else
        print(io, "RawAcqData[$seq_name | $(length(raw.profiles)) Profile(s) of $Nt×$Nc]")
    end
end

"""
Define simulation output for RawData (or RawAcquisitionData)
"""
struct RawData <: SimulationOutput end

function simulation_output(return_type::RawData; kwargs...)
    sim_params_raw = copy(kwargs[:sim_params])
    sim_params_raw["return_type"] = string(kwargs[:sim_params]["return_type"])
    sim_params_raw["sim_method"] = string(kwargs[:sim_params]["sim_method"])
    sim_params_raw["gpu"] = kwargs[:sim_params]["gpu"]
    sim_params_raw["Nthreads"] = kwargs[:sim_params]["Nthreads"]
    sim_params_raw["t_sim_parts"] = kwargs[:t_sim_parts]
    sim_params_raw["type_sim_parts"] = kwargs[:excitation_bool]
    sim_params_raw["Nblocks"] = kwargs[:Nparts]
    sim_params_raw["sim_time_sec"] = kwargs[:timed_tuple_time]
    out = signal_to_raw_data(kwargs[:sig], kwargs[:seq]; phantom_name=kwargs[:obj].name, sys=kwargs[:sys], sim_params=sim_params_raw)
end
