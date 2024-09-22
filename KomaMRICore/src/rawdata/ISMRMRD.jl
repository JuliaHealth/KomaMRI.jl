const b64 = UInt64(1)
const ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1               = 1b64 << ( 01 - 1 )
const ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1                = 1b64 << ( 02 - 1 )
const ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2               = 1b64 << ( 03 - 1 )
const ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2                = 1b64 << ( 04 - 1 )
const ISMRMRD_ACQ_FIRST_IN_AVERAGE                    = 1b64 << ( 05 - 1 )
const ISMRMRD_ACQ_LAST_IN_AVERAGE                     = 1b64 << ( 06 - 1 )
const ISMRMRD_ACQ_FIRST_IN_SLICE                      = 1b64 << ( 07 - 1 )
const ISMRMRD_ACQ_LAST_IN_SLICE                       = 1b64 << ( 08 - 1 )
const ISMRMRD_ACQ_FIRST_IN_CONTRAST                   = 1b64 << ( 09 - 1 )
const ISMRMRD_ACQ_LAST_IN_CONTRAST                    = 1b64 << ( 10 - 1 )
const ISMRMRD_ACQ_FIRST_IN_PHASE                      = 1b64 << ( 11 - 1 )
const ISMRMRD_ACQ_LAST_IN_PHASE                       = 1b64 << ( 12 - 1 )
const ISMRMRD_ACQ_FIRST_IN_REPETITION                 = 1b64 << ( 13 - 1 )
const ISMRMRD_ACQ_LAST_IN_REPETITION                  = 1b64 << ( 14 - 1 )
const ISMRMRD_ACQ_FIRST_IN_SET                        = 1b64 << ( 15 - 1 )
const ISMRMRD_ACQ_LAST_IN_SET                         = 1b64 << ( 16 - 1 )
const ISMRMRD_ACQ_FIRST_IN_SEGMENT                    = 1b64 << ( 17 - 1 )
const ISMRMRD_ACQ_LAST_IN_SEGMENT                     = 1b64 << ( 18 - 1 )
const ISMRMRD_ACQ_IS_NOISE_MEASUREMENT                = 1b64 << ( 19 - 1 )
const ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION             = 1b64 << ( 20 - 1 )
const ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 1b64 << ( 21 - 1 )
const ISMRMRD_ACQ_IS_REVERSE                          = 1b64 << ( 22 - 1 )
const ISMRMRD_ACQ_IS_NAVIGATION_DATA                  = 1b64 << ( 23 - 1 )
const ISMRMRD_ACQ_IS_PHASECORR_DATA                   = 1b64 << ( 24 - 1 )
const ISMRMRD_ACQ_LAST_IN_MEASUREMENT                 = 1b64 << ( 25 - 1 )
const ISMRMRD_ACQ_IS_HPFEEDBACK_DATA                  = 1b64 << ( 26 - 1 )
const ISMRMRD_ACQ_IS_DUMMYSCAN_DATA                   = 1b64 << ( 27 - 1 )
const ISMRMRD_ACQ_IS_RTFEEDBACK_DATA                  = 1b64 << ( 28 - 1 )
const ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 1b64 << ( 29 - 1 )
const ISMRMRD_ACQ_IS_PHASE_STABILIZATION_REFERENCE    = 1b64 << ( 30 - 1 )
const ISMRMRD_ACQ_IS_PHASE_STABILIZATION              = 1b64 << ( 31 - 1 )
const ISMRMRD_ACQ_COMPRESSION1                        = 1b64 << ( 53 - 1 )
const ISMRMRD_ACQ_COMPRESSION2                        = 1b64 << ( 54 - 1 )
const ISMRMRD_ACQ_COMPRESSION3                        = 1b64 << ( 55 - 1 )
const ISMRMRD_ACQ_COMPRESSION4                        = 1b64 << ( 56 - 1 )
const ISMRMRD_ACQ_USER1                               = 1b64 << ( 57 - 1 )
const ISMRMRD_ACQ_USER2                               = 1b64 << ( 58 - 1 )
const ISMRMRD_ACQ_USER3                               = 1b64 << ( 59 - 1 )
const ISMRMRD_ACQ_USER4                               = 1b64 << ( 60 - 1 )
const ISMRMRD_ACQ_USER5                               = 1b64 << ( 61 - 1 )
const ISMRMRD_ACQ_USER6                               = 1b64 << ( 62 - 1 )
const ISMRMRD_ACQ_USER7                               = 1b64 << ( 63 - 1 )
const ISMRMRD_ACQ_USER8                               = 1b64 << ( 64 - 1 )
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
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}(), ndims=2
)
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
    #userParameters <- sim_params
    for (key, val) in sim_params
        if typeof(val) <: Integer #Fixes problem with bools
            sim_params[key] = Int(val)
        end
    end
    #XML header
    params = Dict(
        #AcquisitionSystemInformation
        "systemVendor"                   => "KomaMRI.jl", #String
        "systemModel"                    => string(pkgversion(@__MODULE__)), #String
        "systemFieldStrength_T"          => sys.B0, #Float
        "institutionName"                => "Pontificia Universidad Catolica de Chile", #String
        #subjectInformation
        "patientName"                    => phantom_name,
        #experimentalConditions
        "H1resonanceFrequency_Hz"        => floor(Int64, γ * sys.B0), #Long (Int)
        #measurementInformation
        "protocolName"                   => haskey(seq.DEF,"Name") ? seq.DEF["Name"] : "NoName", #String
        # "trajectoryDescription"          => Dict{String, Any}("comment"=>""), #You can put wathever you want here: comment, bandwidth, MaxGradient_G_per_cm, MaxSlewRate_G_per_cm_per_s, interleaves, etc
        #encoding
        #   encodedSpace
        "encodedSize"                    => [Nx, Ny, 1],                        #encodedSpace>matrixSize
        "encodedFOV"                     => Float32.([FOVx, FOVy, 1e-3]*1e3),   #encodedSpace>fieldOfView_mm
        #   reconSpace
        "reconSize"                      => [Nx+Nx%2, Ny+Ny%2, 1],              #reconSpace>matrixSize
        "reconFOV"                       => Float32.([FOVx, FOVy, 1e-3]*1e3),   #reconSpace>fieldOfView_mm
        #encodingLimits
        "enc_lim_kspace_encoding_step_0" => Limit(0, Nx-1, ceil(Int, Nx / 2)),  #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_1" => Limit(0, Ny-1, ceil(Int, Ny / 2)),  #min, max, center, e.g. partition encoding number
        "enc_lim_kspace_encoding_step_2" => Limit(0, 0, 0),                     #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),                     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, 0, 0),                     #min, max, center, e.g. imaging slice number
        "enc_lim_contrast"               => Limit(0, 0, 0),                     #min, max, center, e.g. echo number in multi-echo
        "enc_lim_phase"                  => Limit(0, 0, 0),                     #min, max, center, e.g. cardiac phase number
        "enc_lim_repetition"             => Limit(0, 0, 0),                     #min, max, center, e.g. dynamic number for dynamic scanning
        "enc_lim_set"                    => Limit(0, 0, 0),                     #min, max, center, e.g. flow encoding set
        "enc_lim_segment"                => Limit(0, 0, 0),                     #min, max, center, e.g. segment number for segmented acquisition
        "trajectory"                     => "other",
        #sequenceParameters
        # "TR"                             => 0,
        # "TE"                             => 0,
        # "TI"                             => 0,
        # "flipAngle_deg"                  => 0,
        # "echo_spacing"                   => 0,
        "userParameters"                 => sim_params, #Dict with parameters
    )

    #Then, we define the Profiles
    profiles = Profile[]
    t_acq = get_adc_sampling_times(seq)
    Nadcs = sum(is_ADC_on.(seq))
    NadcsPerImage = floor(Int, Nadcs / Nz)
    scan_counter = 0
    nz = 0
    current = 1
    for s = seq #Iterate over sequence blocks
        if is_ADC_on(s)
            Nsamples = s.ADC.N[1]
            Δt_us = floor( s.ADC.T[1] / (Nsamples - 1) * 1e6 )
            t0_us = floor( t_acq[current] * 1e6 )
            flag  = 0
            if scan_counter == 0
                flag += ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1
                flag += ISMRMRD_ACQ_FIRST_IN_SLICE
            elseif scan_counter == Nadcs - 1
                flag += ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1
                flag += ISMRMRD_ACQ_LAST_IN_SLICE
            end
            #Header of profile data, head::AcquisitionHeader
            head = AcquisitionHeader(
                UInt16(1), #version uint16: First unsigned int indicates the version
                UInt64(flag), #flags uint64: bit field with flags, noise measure, calibration, coil sens, measure, etc
                UInt32(0), #measurement_uid uint32: Unique ID for the measurement
                UInt32(scan_counter), #scan_counter uint32: Current acquisition number in the measurement
                UInt32(t0_us), #acquisition_time_stamp uint32: Acquisition clock, I am "miss"-using this variable to store t0 in us
                UInt32.((0, 0, 0)), #physiology_time_stamp uint32x3: Physiology time stamps, e.g. ecg, breating, etc.
                UInt16(Nsamples), #number_of_samples uint16
                UInt16(1), #available_channels uint16: Available coils
                UInt16(1), #active_channels uint16: Active coils on current acquisiton
                Tuple(UInt64(0) for i=1:16), #channel_mask uint64x16: Active coils on current acquisiton
                UInt16(0), #discard_pre uint16: Samples to be discarded at the beginning of acquisition
                UInt16(0), #discard_post uint16: Samples to be discarded at the end of acquisition
                UInt16(0), #center_sample uint16: Sample at the center of k-space
                UInt16(0), #encoding_space_ref uint16: Reference to an encoding space, typically only one per acquisition
                UInt16(ndims), #trajectory_dimensions uint16: Indicates the dimensionality of the trajectory vector (0 means no trajectory)
                Float32(Δt_us), #sample_time_us float32: Time between samples in micro seconds, sampling BW
                Float32.((0, 0, 0)), #position float32x3: Three-dimensional spatial offsets from isocenter
                Float32.((1, 0, 0)), #read_dir float32x3: Directional cosines of the readout/frequency encoding
                Float32.((0, 1, 0)), #phase_dir float32x3: Directional cosines of the phase
                Float32.((0, 0, 1)), #slice_dir float32x3: Directional cosines of the slice direction
                Float32.((0, 0, 0)), #patient_table_position float32x3: Patient table off-center
                EncodingCounters( #idx uint16x17: Encoding loop counters
                    UInt16(scan_counter), #kspace_encode_step_1 uint16: e.g. phase encoding line number
                    UInt16(0), #kspace_encode_step_2 uint16: e.g. partition encoding number
                    UInt16(0), #average uint16: e.g. signal average number
                    UInt16(nz), #slice uint16: e.g. imaging slice number
                    UInt16(0), #contrast uint16: e.g. echo number in multi-echo
                    UInt16(0), #phase uint16: e.g. cardiac phase number
                    UInt16(0), #repetition uint16: e.g. dynamic number for dynamic scanning
                    UInt16(0), #set uint16: e.g. flow encoding set
                    UInt16(0), #segment uint16: e.g. segment number for segmented acquisition
                    Tuple(UInt16(0) for i=1:8) #user uint16x8: Free user parameters
                ),
                Tuple(Int32(0) for i=1:8), #user_int int32x8: Free user parameters
                Tuple(Float32(0) for i=1:8) #user_float float32x8: Free user parameters
            )
            #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples
            traj = ktraj[1:ndims, current:current+Nsamples-1]
            #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
            dat =  signal[current:current+Nsamples-1, :]
            #Saving profile
            push!(profiles, Profile(head, Float32.(traj), ComplexF32.(dat)))
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

Base.:+(sig1::RawAcquisitionData, sig2::RawAcquisitionData) = RawAcquisitionData(
    sig1.params,
    [Profile(
        sig1.profiles[i].head,
        sig1.profiles[i].traj,
        sig1.profiles[i].data .+ sig2.profiles[i].data
    ) for i=1:length(sig1.profiles)]
)

Base.isapprox(sig1::RawAcquisitionData, sig2::RawAcquisitionData; kwargs...) = begin
    (length(sig1.profiles) == length(sig2.profiles)) || return false

    for i=1:length(sig1.profiles)
        isapprox(sig1.profiles[i].data, sig2.profiles[i].data; kwargs...) || return false
    end

    return true
end
