ifftc(x;dims=[1,2])=fftshift(ifft(ifftshift(x,dims),dims),dims)*prod(size(x)[dims])
fftc(x;dims=[1,2]) =fftshift(fft(ifftshift(x,dims),dims),dims)/prod(size(x, d) for d in dims)

# Allowed actions, default first. See DEFAULT_RECON_POLICY below for their meanings.
const _RECON_POLICY_OPTIONS = (
    ROLE=(:separate,),
    DIM=(:separate,),
    LIN=(:encoding,),
    PAR=(:encoding,),
    SLC=(:separate,),
    SEG=(:combine, :separate),
    REP=(:separate,),
    AVG=(:separate, :combine),
    SET=(:separate,),
    ECO=(:axis, :separate),
    PHS=(:separate,),
    ACQ=(:ignore,),
    TRID=(:ignore,),
    COIL=(:axis, :rss),
)

"""
    DEFAULT_RECON_POLICY

Default label handling for [`reconstruct_with_labels`](@ref). Pass a partial named tuple
such as `recon_policy=(AVG=:combine, COIL=:rss)` to override the configurable actions.

| Dimension | Meaning | Default | Alternative |
|:----------|:--------|:--------|:------------|
| `ROLE` | Acquisition purpose, derived from MRD flags | `:separate` | — |
| `DIM` | Stored trajectory coordinate count, not image rank | `:separate` | — |
| `LIN` | Phase-encoding line index | `:encoding` | — |
| `PAR` | Partition-encoding index within a volume | `:encoding` | — |
| `SLC` | Separate acquired slice | `:separate` | — |
| `SEG` | Segments contributing to an image | `:combine` | `:separate` |
| `REP` | Repetition counter | `:separate` | — |
| `AVG` | Repeated acquisition for averaging | `:separate` | `:combine` |
| `SET` | Acquisition set counter | `:separate` | — |
| `ECO` | Echo label, stored as MRD contrast | `:axis` | `:separate` |
| `PHS` | Physiological phase counter, not complex signal phase | `:separate` | — |
| `ACQ` | Acquisition counter, stored in MRD `idx.user[1]` | `:ignore` | — |
| `TRID` | Interpreter-specific counter, in MRD `idx.user[2]` | `:ignore` | — |
| `COIL` | Receiver channels, not an ADC label | `:axis` | `:rss` |

Actions describe both the MRIReco call boundary and the returned image organization:

- `:separate`: split profiles by every observed combination of the separate labels.
  Each combination gets its own MRIReco call and entry in `result.images`; absent
  combinations are not created.
- `:axis`: keep the values together in one MRIReco call and retain its echo or coil
  axis in one image `AxisArray`. Values are not averaged or otherwise combined.
- `:encoding`: retain `LIN`/`PAR` as k-space encoding indices within an image, not as
  separate images. These indices do not determine stored trajectory dimensionality.
- `:combine`: `AVG` sums the complex signals of matching profiles with identical
  trajectories; `SEG` joins segment profiles for one reconstruction without averaging.
- `:ignore`: do not split or combine using this counter. Original counter values are
  still retained in each entry's `source` metadata.
- `:rss`: replace reconstructed coil images with `sqrt(sum(abs2, coils))`. This is
  magnitude-only combination, without coil-sensitivity estimation or retained phase.

`ROLE` is selected by MRD flag precedence: noise, phase correction, navigator,
calibration-and-imaging, calibration, phase-stabilization reference, phase stabilization,
feedback, dummy, then surface-coil correction. Unflagged and calibration-and-imaging
profiles are imaging. Imaging, navigator, phase-correction, calibration, and both
phase-stabilization roles are reconstructable; the other roles are excluded.
"""
const DEFAULT_RECON_POLICY = map(first, _RECON_POLICY_OPTIONS)
const ReconstructionResult = NamedTuple{(:policy, :images)}

"""
    reconstruction_policy(overrides=NamedTuple())

Validate and merge the partial `overrides` named tuple with [`DEFAULT_RECON_POLICY`](@ref).
See that policy's documentation for every dimension, allowed action, and default.
"""
function reconstruction_policy(overrides=NamedTuple())
    for (name, action) in pairs(overrides)
        action in getproperty(_RECON_POLICY_OPTIONS, name) || throw(ArgumentError(
            "Unsupported $name reconstruction policy $action; expected one of " *
            join(getproperty(_RECON_POLICY_OPTIONS, name), ", ")
        ))
    end
    return merge(DEFAULT_RECON_POLICY, overrides)
end

# ISMRMRD flags are authoritative. The first matching flag sets the exclusive role;
# calibration-and-imaging data is explicitly treated as imaging.
const _ACQUISITION_ROLE_FLAGS = (
    :noise => KomaMRICore.ISMRMRD_ACQ_IS_NOISE_MEASUREMENT,
    :phase_correction => KomaMRICore.ISMRMRD_ACQ_IS_PHASECORR_DATA,
    :navigator => KomaMRICore.ISMRMRD_ACQ_IS_NAVIGATION_DATA,
    :imaging => KomaMRICore.ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING,
    :calibration => KomaMRICore.ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION,
    :phase_stabilization_reference =>
        KomaMRICore.ISMRMRD_ACQ_IS_PHASE_STABILIZATION_REFERENCE,
    :phase_stabilization => KomaMRICore.ISMRMRD_ACQ_IS_PHASE_STABILIZATION,
    :feedback => KomaMRICore.ISMRMRD_ACQ_IS_HPFEEDBACK_DATA |
        KomaMRICore.ISMRMRD_ACQ_IS_RTFEEDBACK_DATA,
    :dummy => KomaMRICore.ISMRMRD_ACQ_IS_DUMMYSCAN_DATA,
    :surface_coil => KomaMRICore.ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA,
)

function _acquisition_role(flags)
    for (role, flag) in _ACQUISITION_ROLE_FLAGS
        !iszero(flags & flag) && return role
    end
    return :imaging
end

_reconstructable(role) = role in (
    :imaging,
    :navigator,
    :phase_correction,
    :calibration,
    :phase_stabilization_reference,
    :phase_stabilization,
)

function _profile_dimensions(profile)
    idx = profile.head.idx
    return (
        ROLE=_acquisition_role(profile.head.flags),
        DIM=profile.head.trajectory_dimensions,
        LIN=idx.kspace_encode_step_1,
        PAR=idx.kspace_encode_step_2,
        SLC=idx.slice,
        SEG=idx.segment,
        REP=idx.repetition,
        AVG=idx.average,
        SET=idx.set,
        ECO=idx.contrast,
        PHS=idx.phase,
        ACQ=idx.user[1],
        TRID=idx.user[2],
    )
end

function _batch_key(profile, names)
    dimensions = _profile_dimensions(profile)
    return NamedTuple{names}(Tuple(getproperty(dimensions, name) for name in names))
end

const _MRD_COUNTERS = (
    SLC=:slice,
    SEG=:segment,
    REP=:repetition,
    AVG=:average,
    SET=:set,
    ECO=:contrast,
    PHS=:phase,
)

function _prepare_batch_profile(profile, policy)
    head = deepcopy(profile.head)
    if !isempty(profile.traj)
        head.flags &= ~KomaMRICore.ISMRMRD_ACQ_IS_REVERSE
    end
    for (name, field) in pairs(_MRD_COUNTERS)
        getproperty(policy, name) === :separate && setproperty!(head.idx, field, 0)
    end
    return Profile(head, profile.traj, profile.data)
end

function _average_profile_key(profile)
    idx = profile.head.idx
    is_reverse = !iszero(profile.head.flags & KomaMRICore.ISMRMRD_ACQ_IS_REVERSE)
    return (
        profile.head.trajectory_dimensions,
        idx.kspace_encode_step_1,
        idx.kspace_encode_step_2,
        idx.slice,
        idx.contrast,
        idx.phase,
        idx.repetition,
        idx.set,
        idx.segment,
        is_reverse,
    )
end

function _combine_averages(profiles)
    averages = sort!(unique(profile.head.idx.average for profile in profiles))
    length(averages) == 1 && return profiles

    groups = [
        filter(profile -> profile.head.idx.average == average, profiles) for
        average in averages
    ]
    profile_keys = _average_profile_key.(first(groups))
    all(group -> _average_profile_key.(group) == profile_keys, groups) ||
        error("Cannot combine averages with different acquisition profiles")

    return map(eachindex(profile_keys)) do index
        group = getindex.(groups, index)
        reference = first(group)
        all(profile -> profile.traj == reference.traj, group) ||
            error("Cannot combine averages with different trajectories")
        head = deepcopy(reference.head)
        head.idx.average = 0
        data = reduce(+, getproperty.(group, :data))
        Profile(head, reference.traj, data)
    end
end

_retained_samples(profile) =
    (1 + profile.head.discard_pre):(size(profile.data, 1) - profile.head.discard_post)

function _encoding_axes(raw)
    samples = sum(profile -> length(_retained_samples(profile)), raw.profiles)
    normalization = get(get(raw.params, "userParameters", Dict()), "KomaTrajectoryScale", nothing)
    # Exported data has a physical tolerance; foreign MRD falls back to Float32 roundoff.
    tolerance = isnothing(normalization) ? eps(Float32) * _trajectory_scale(raw) :
        KomaMRICore.KSPACE_DIMENSION_TOLERANCE / normalization
    return filter(1:first(raw.profiles).head.trajectory_dimensions) do dimension
        lower = minimum(profile -> minimum(@view(profile.traj[dimension, _retained_samples(profile)])), raw.profiles)
        upper = maximum(profile -> maximum(@view(profile.traj[dimension, _retained_samples(profile)])), raw.profiles)
        mean_coordinate = sum(raw.profiles) do profile
            sum(Float64, @view(profile.traj[dimension, _retained_samples(profile)]))
        end / samples
        max(upper - mean_coordinate, mean_coordinate - lower) > tolerance
    end
end

function _reconstruction_geometry(raw)
    dimensions = first(raw.profiles).head.trajectory_dimensions
    dimensions == 0 && return raw
    encoding_axes = _encoding_axes(raw)
    isempty(encoding_axes) && return raw
    if length(encoding_axes) > 1 && any(d -> raw.params["encodedSize"][d] <= 1, encoding_axes)
        throw(ArgumentError(
            "Cannot determine the reconstruction size. " *
            "Add Nx/Ny/Nz to the .seq file's [DEFINITIONS] section."
        ))
    end
    length(encoding_axes) == dimensions && dimensions > 1 && return raw

    # Geometry belongs to this role/policy batch, not to the complete MRD acquisition.
    params = copy(raw.params)
    order = [encoding_axes; setdiff(1:3, encoding_axes)]
    for name in ("encodedSize", "reconSize", "encodedFOV", "reconFOV")
        haskey(params, name) && (params[name] = raw.params[name][order])
    end
    for name in ("encodedSize", "reconSize")
        params[name][length(encoding_axes)+1:3] .= 1
    end
    if length(encoding_axes) == 1
        # A navigator inherits the imaging matrix in a mixed acquisition. Its readout
        # sample count supplies the line size; standalone 1D data keeps its own matrix.
        if count(>(1), raw.params["encodedSize"]) > 1 || params["encodedSize"][1] == 1
            samples = maximum(profile -> length(_retained_samples(profile)), raw.profiles)
            params["encodedSize"][1] = samples
            params["reconSize"][1] = samples + samples % 2
        end
        params["reconSize"][2] = 2
    end
    profiles = map(raw.profiles) do profile
        head = deepcopy(profile.head)
        head.trajectory_dimensions = max(length(encoding_axes), 2)
        head.idx.kspace_encode_step_2 = 0
        trajectory = profile.traj[encoding_axes, :]
        if length(encoding_axes) == 1
            head.idx.kspace_encode_step_1 = 0
            trajectory = vcat(trajectory, zeros(Float32, 1, size(trajectory, 2)))
        end
        Profile(head, trajectory, profile.data)
    end
    return RawAcquisitionData(params, profiles)
end

function _source_labels(profiles)
    names = (:LIN, :PAR, :SLC, :SEG, :REP, :AVG, :SET, :ECO, :PHS, :ACQ, :TRID)
    dimensions = _profile_dimensions.(profiles)
    values = map(names) do name
        sort!(unique(getproperty.(dimensions, name)))
    end
    return NamedTuple{names}(values)
end

function reconstruction_batches(raw, policy)
    isempty(raw.profiles) && return NamedTuple[]
    separate_labels = Tuple(name for (name, action) in pairs(policy) if action === :separate)
    label_type = typeof(_batch_key(first(raw.profiles), separate_labels))
    grouped_profiles = Dict{label_type,Vector{Profile}}()
    for profile in raw.profiles
        labels = _batch_key(profile, separate_labels)
        _reconstructable(labels.ROLE) || continue
        push!(get!(Vector{Profile}, grouped_profiles, labels), profile)
    end
    batch_keys = sort!(collect(keys(grouped_profiles)); by=labels ->
        (labels.ROLE === :imaging ? 0 : 1, Tuple(labels)...))
    return [
        begin
            profiles = grouped_profiles[batch_key]
            source = _source_labels(profiles)
            batch_profiles = policy.AVG === :combine ?
                _combine_averages(profiles) : profiles
            batch = RawAcquisitionData(
                copy(raw.params),
                [_prepare_batch_profile(profile, policy) for profile in batch_profiles],
            )
            (; labels=batch_key, source, raw=_reconstruction_geometry(batch))
        end for batch_key in batch_keys
    ]
end

function _trajectory_scale(raw)
    scale = mapreduce(
        profile -> isempty(profile.traj) ? 0.0f0 :
            2 * maximum(abs, @view(profile.traj[:, _retained_samples(profile)])),
        max,
        raw.profiles;
        init=0.0f0,
    )
    return iszero(scale) ? one(scale) : scale
end

function _prepare_acquisition!(acquisition, scale)
    for trajectory in acquisition.traj
        trajectory.circular = false
        trajectory.nodes ./= scale
    end
    return nothing
end

function _rss(image)
    data = sqrt.(sum(abs2, image; dims=5))
    image_axes = AxisArrays.axes(image)
    return AxisArrays.AxisArray(
        data,
        image_axes[1:4]...,
        AxisArrays.Axis{:coils}((:rss,)),
        image_axes[6],
    )
end

function _reconstruct_batch(batch, rec_params, coil_policy)
    raw = batch.raw
    acquisition = AcquisitionData(raw)
    _prepare_acquisition!(acquisition, _trajectory_scale(raw))
    dimensions = size(first(acquisition.traj).nodes, 1)
    params = Dict{Symbol,Any}(pairs(rec_params))
    params[:reconSize] = Tuple(Int.(raw.params["reconSize"][1:dimensions]))
    get!(params, :densityWeighting, true)
    image = reconstruction(acquisition, params)
    # Do not expose the padded phase-encoding row of a 1D acquisition.
    dimensions == 2 && raw.params["encodedSize"][2] == 1 &&
        (image = image[:, 1:1, :, :, :, :])
    coil_policy === :rss && (image = _rss(image))
    return (; batch.labels, batch.source, image)
end

_magnitude_group(entry, echo) = (
    entry.labels.ROLE, entry.labels.DIM, get(entry.labels, :SET, nothing), echo,
)

function _with_magnitude_limits(images)
    first_group = _magnitude_group(first(images), first(first(images).source.ECO))
    limits = Dict{typeof(first_group),Tuple{Float64,Float64}}()
    for entry in images, (contrast, echo) in enumerate(entry.source.ECO)
        group = _magnitude_group(entry, echo)
        zmin, zmax = extrema(abs, @view(entry.image[:, :, :, contrast, :, :]))
        scale = prod(size(entry.image)[1:2])
        lower, upper = get(limits, group, (Inf, -Inf))
        limits[group] = (min(lower, zmin * scale), max(upper, zmax * scale))
    end
    return [
        (; entry..., magnitude_limits=[limits[_magnitude_group(entry, echo)] for echo in entry.source.ECO])
        for entry in images
    ]
end

"""
    reconstruct_with_labels(
        raw,
        ;
        recon_policy=NamedTuple(),
        rec_params=Dict{Symbol,Any}(),
    )

Reconstruct MRD data with MRIReco after applying a policy to its acquisition labels.
`recon_policy` contains optional overrides of [`DEFAULT_RECON_POLICY`](@ref), whose
documentation defines every dimension, action (`:separate`, `:axis`, `:encoding`,
`:combine`, `:ignore`, `:rss`), and acquisition role. `rec_params` supplies MRIReco's
reconstruction options; each batch's MRD `reconSize` determines its image matrix.

Profiles are separated by role and policy before determining geometry or normalizing
trajectories. Only coordinate variation within that batch determines its encoding axes;
storing three k-space coordinates does not make every role a 3D acquisition. Axis-aligned
1D readouts in x, y, or z use MRIReco's padded 2D path internally and return a single line.
For a line within a multi-dimensional acquisition, the retained ADC sample count supplies
its matrix size instead of the imaging matrix. Geometry and scaling exclude MRD samples
marked for discarding, as MRIReco does. The input MRD data is not modified.

With other labels fixed, `SLC=:separate` and `REP=:separate` turn observed `(SLC, REP)`
pairs `(0, 0)`, `(0, 1)`, and `(1, 0)` into three independent images, not four. KomaUI
provides SLC and REP controls and only offers combinations present in the data.

With all other labels fixed and two ECO values, `ECO=:axis` produces one `AxisArray`
with two echoes through one MRIReco call; `ECO=:separate` produces two entries through
two calls. Both provide an ECO control in KomaUI. Neither action averages the echoes.
Likewise, `AVG=:separate` keeps individual acquisitions, while `AVG=:combine` sums
their matching complex data before reconstruction; `COIL=:rss` instead combines coil
magnitudes after reconstruction.

The result is a named tuple `(policy, images)`. Each element of `images` is a named tuple
`(labels, source, image, magnitude_limits)`, where `image` is MRIReco's `AxisArray`.
`magnitude_limits` contains one display-range pair per `source.ECO` value, computed once
after reconstruction. Ranges are shared across slices, partitions, averages, repetitions,
phases, segments, and coils, but separate for each `(ROLE, DIM, SET, ECO)` combination.
They use KomaUI's magnitude scaling by the first two image dimensions; phase and k-space
display ranges remain fixed. Treat the returned images as immutable or recompute their
limits after changing image data.
"""
function reconstruct_with_labels(
    raw;
    recon_policy=NamedTuple(),
    rec_params=Dict{Symbol,Any}(),
)
    policy = reconstruction_policy(recon_policy)
    batches = reconstruction_batches(raw, policy)
    isempty(batches) && error("Raw data has no reconstructable profiles")
    images = _with_magnitude_limits([
        _reconstruct_batch(batch, rec_params, policy.COIL) for batch in batches
    ])
    return (; policy, images)
end
