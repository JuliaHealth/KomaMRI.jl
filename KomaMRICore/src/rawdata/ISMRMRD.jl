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
    raw = signal_to_raw_data(signal, seq; phantom_name, sys, sim_params, ndims=2, use_ndseq=false)

Transforms the raw signal into a RawAcquisitionData struct (nearly equivalent to the ISMRMRD
format) used for reconstruction with MRIReco.

# Arguments
- `signal`: (`::Matrix{Complex}`) raw signal matrix
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `phantom_name`: (`::String`, `="Phantom"`) phantom name
- `sys`: (`::Scanner`, `=Scanner()`) Scanner struct
- `sim_params`: (`::Dict{String, Any}`, `=Dict{String,Any}()`) simulation parameter dictionary
- `ndims` : (`::Integer`, `=2`) number of dimensions of the reconstruction
- `use_ndseq` = (`Bool`, `=false`) attempts to estimate dimension of reconstruction from seq file

# Returns
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_se.seq")

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> sim_params = KomaMRICore.default_sim_params(); sim_params["return_type"] = "mat"

julia> signal = simulate(obj, seq, sys; sim_params)

julia> raw = signal_to_raw_data(signal, seq; ndims=3)

julia> plot_signal(raw)
```
"""
function signal_to_raw_data(
    signal, seq;
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}(), ndims=2, use_ndseq=true
)

    #Number of samples and FOV
    Nd_seq, Nx, Ny, Nz, Ns, FOVx, FOVy, FOVz, Δx, ktraj = estimate_seq_recon_dimension(seq; sim_params)
    if use_ndseq
        if ndims != Nd_seq; @warn("Using estimated dimension of $Nd_seq for .mrd file."); ndims = Nd_seq; end
    else
        if ndims != Nd_seq; @warn("Seqfile is $Nd_seq dimensional but .mrd file will be $ndims."); end
    end
    if ndims == 2
        @info "Creating 2D ISMRMRD file ..."
    elseif ndims == 3        
        @info "Creating 3D ISMRMRD file ..."
    else
        @info("Creating a $ndims dimensional ISMRMRD file...")
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
        "encodedSize"                    => [Nx, Ny, Nz],                       #encodedSpace>matrixSize
        "encodedFOV"                     => Float32.([FOVx, FOVy, FOVz]*1e3),   #encodedSpace>fieldOfView_mm
        #   reconSpace
        "reconSize"                      => [Nx, Ny, Nz],                       #reconSpace>matrixSize
        "reconFOV"                       => Float32.([FOVx, FOVy, FOVz]*1e3),   #reconSpace>fieldOfView_mm
        #encodingLimits
        "enc_lim_kspace_encoding_step_0" => Limit(0, Nx-1, floor(Int, Nx / 2)), #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_1" => Limit(0, Ny-1, floor(Int, Ny / 2)), #min, max, center, e.g. partition encoding number
        "enc_lim_kspace_encoding_step_2" => Limit(0, Nz-1, floor(Int, Nz / 2)), #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),                     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, Ns-1, floor(Int, Ns / 2)), #min, max, center, e.g. imaging slice number
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
    Nro = sum(is_ADC_on.(seq))
    #Nyt = max(Ny, 1); Nst = max(Ns, 1); Nzt = max(Nz, 1) #zero is really one...move to estimate_seq_recon_dimension
    NroPerPE1   = ceil(Int, Nro / (Ny*Nz*Ns))
    NroPerSlice = ceil(Int, Nro / Ns)
    NroPerPE2   = ceil(Int, Nro / Nz)
    @debug "Nro, NroPerPE1, NroPerSlice, NroPerPE2 = $Nro, $NroPerPE1, $NroPerSlice, $NroPerPE2"
    scan_counter = 0
    ny = 0 #PE1 counter
    ns = 0 #Slice counter
    nz = 0 #PE2 counter
    current = 1 #used for traj and data partition
    for s = seq #Iterate over sequence blocks
        if is_ADC_on(s)
            Nsamples = s.ADC.N[1]
            Δt_us = floor( s.ADC.T[1] / (Nsamples - 1) * 1e6 )
            t0_us = floor( t_acq[current] * 1e6 )
            flag  = 0
            if scan_counter == 0
                flag += ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1
                flag += ISMRMRD_ACQ_FIRST_IN_SLICE
                if Nz > 1; flag += ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2; end
            elseif scan_counter == Nro - 1 #Needs fixing? CAC 240609 ***
                flag += ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1
                flag += ISMRMRD_ACQ_LAST_IN_SLICE
                if Nz > 1; flag += ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2; end
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
                    UInt16(ny), #kspace_encode_step_1 uint16: e.g. phase encoding line number
                    UInt16(nz), #kspace_encode_step_2 uint16: e.g. partition encoding number
                    UInt16(0),  #average uint16: e.g. signal average number
                    UInt16(ns), #slice uint16: e.g. imaging slice number
                    UInt16(0),  #contrast uint16: e.g. echo number in multi-echo
                    UInt16(0),  #phase uint16: e.g. cardiac phase number
                    UInt16(0),  #repetition uint16: e.g. dynamic number for dynamic scanning
                    UInt16(0),  #set uint16: e.g. flow encoding set
                    UInt16(0),  #segment uint16: e.g. segment number for segmented acquisition
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
            scan_counter += 1 #another ro
            current += Nsamples
            if scan_counter % NroPerPE1 == 0
                ny += 1 #another kspace_encode_step_1
            end
            ny = ny % Ny
            if scan_counter % NroPerPE2 == 0
                nz += 1 #another image
            end
            nz = nz % Nz
            if  scan_counter % NroPerSlice == 0
                ns += 1 #another slice
            end
            ns = ns % Ns
            # more here for other counters, i.e. dynamic
        end
    end
    @debug "scan_counter, ny, nz, ns = $scan_counter, $ny, $nz, $ns"
    return RawAcquisitionData(params, profiles)
end

"""
    Nd_seq, Nx, Ny, Nz, Ns, FOVx, FOVy, FOVz, Δx, ktraj = estimate_seq_recon_dimension(seq; sim_params)

Utility function for the best estimate of the reconstruction dimension.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `sim_params`: (`::Dict{String, Any}`, `=Dict{String,Any}()`) simulation parameter dictionary (IGNORED)

# Returns
- `Nd_seq`: (`Int`) Estimated reconstruction dimension of seq.

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_se.seq")

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> sim_params = KomaMRICore.default_sim_params(); sim_params["return_type"] = "mat"

julia> Nd_seq = estimate_seq_recon_dimension(seq; sim_params)
```
"""
function estimate_seq_recon_dimension(seq; sim_params=Dict{String,Any}(),
)
    #remove signal and sim_params...not needed   
    
    #K-space info
    _, ktraj = get_kspace(seq) #kspace information
    mink = minimum(ktraj, dims=1)
    maxk = maximum(ktraj, dims=1)
    Wk = maxk .- mink
    idxs_zero = findall(iszero, Wk)  #check for zeros
    Wk_el_iszero = iszero.( Wk) # bool array for later
    Wk[idxs_zero] .= 1.0e-6
    Δx = 1 ./ Wk[1:3] #[m] x-y-z
    k_radius= mapslices(norm, ktraj, dims=2) #limit to first N readouts?
    k_radius_norm = k_radius ./ maximum(k_radius)

    #estimate sequence acquisition structure
    Ns_seq = length(unique(seq.RF.Δf)) #slices, slabs, or Cartesean off-center, preps need to & with ADC *** CAC 240708
    Np_seq = maximum(adc.N for adc = seq.ADC) #readout length
    Nv_seq = sum(map(is_ADC_on, seq)) #total number of readouts or views
    @debug "Ns_seq, Np_seq, Nv_seq = $Ns_seq, $Np_seq, $Nv_seq"
    
    #estimate FOV and N from sequence and K-space info
    FOVx_k = Np_seq*Δx[1]
    FOVy_k = Np_seq*Δx[2]
    FOVz_k = Np_seq*Δx[3]
    Nx_k = ceil(Int64, FOVx_k / Δx[1])
    Ny_k = ceil(Int64, FOVy_k / Δx[2])
    Nz_k = ceil(Int64, FOVz_k / Δx[3])
    @debug "FOVx_k, FOVy_k, FOVz_k = $FOVx_k, $FOVy_k, $FOVz_k"
    @debug "Nx_k, Ny_k, Nz_k = $Nx_k, $Ny_k, $Nz_k"
    
    # some guesses
    seq_2d=false; seq_3d=false; seq_cartesean=false; seq_radial=false; seq_spiral=false
    if FOVz_k > 100
        seq_2d=true
    elseif FOVz_k > 0
        seq_3d=true
    else
        @warn "This seq file does not seem to be an imaging sequence."    
    end
    if Np_seq > 5*Nv_seq
        seq_spiral = true
    elseif sum( k_radius_norm .< 3/Np_seq) >= Nv_seq #Center of K-space highly sampled
        seq_radial=true
    else
        seq_cartesean=true
    end
    @debug "seq_2d, seq_3d, seq_cartesean, seq_radial, seq_spiral = $seq_2d, $seq_3d, $seq_cartesean, $seq_radial, $seq_spiral"
    
    # Guessing Cartesean recon dimensions
    # ideally all estimates of recon dimensions in one place, as late as possible, non-Cartesean ?? *** CAC 240708
    Ns = Int64(get(seq.DEF, "Ns", Ns_seq)) #slices or slabs
    Nx = Int64(get(seq.DEF, "Nx", Np_seq))
    if seq_cartesean
        if seq_2d
            Ny = Int64(get(seq.DEF, "Ny", ceil(Int64, Nv_seq/Ns))) #pe1
            Nz = 1
        elseif seq_3d
            Ny = Int64(get(seq.DEF, "Ny", Ny_k)) #pe1
            Nz = Int64(get(seq.DEF, "Nz", ceil(Int64, Nv_seq/Ny_k))) #pe2
        end
    else
        @warn "Non-Cartesean recon is still under development."
        Ny = ceil(Int64, Nv_seq)
        if Ny == 1 Ns = 1 end #sinle shot spiral
        Nz = 1
    end
    
    # explicit reads from seq.DEF
    #Nx = ceil(Int64, get(seq.DEF, "Nx", 1))
    #Ny = ceil(Int64, get(seq.DEF, "Ny", 1))
    #Nz = ceil(Int64, get(seq.DEF, "Nz", 1))
    #Ns = ceil(Int64, get(seq.DEF, "Ns", 1)) # number of slices or slabs
    
    #if Nx == 1  Nx = ceil(Int64, FOVx / Δx[1])  end
    #if Ny == 1  Ny = ceil(Int64, FOVy / Δx[2])  end
    #if Nz == 1  Nz = ceil(Int64, FOVz / Δx[3])  end  

    if haskey(seq.DEF, "FOV")   #needs info or warning for other substitutions???
        FOVx, FOVy, FOVz = seq.DEF["FOV"] #[m]
        if FOVx > 1.0 #do this in steps until < 1?? seems like mm or cm possible
            FOVx *= 1e-2
            FOVy *= 1e-2
            FOVz *= 1e-2
            @warn "Scaling FOV to m from older pulseq file cm."
        end
        if seq_3d
            FOVz = FOVz
        else
            FOVz = 0.0
        end  
    elseif seq_cartesean
        @warn "Estimating FOV parameters from k-space."
        FOVx = FOVx_k
        FOVy = FOVy_k
        if seq_3d
            FOVz = FOVz_k
        else
            FOVz = 0.0
        end
    end
    Nd_seq = (FOVx > .05) + (FOVy > .05) + (FOVz > .05)
    s_ktraj = size( ktraj)
    @debug "Nd_seq = $Nd_seq"
    @debug "Nx, Ny, Nz, Ns = $Nx, $Ny, $Nz, $Ns"
    @debug "FOVx, FOVy, FOVz = $FOVx, $FOVy, $FOVz"
    @debug "Δx, size( ktraj) = $Δx, $s_ktraj"
    return Nd_seq, Nx, Ny, Nz, Ns, FOVx, FOVy, FOVz, Δx, ktraj
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
