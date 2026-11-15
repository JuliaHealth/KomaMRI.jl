function combine_coil_signal!(sig, Mxy::AbstractMatrix, sens, ::KA.CPU)
    sig .= transpose(Mxy) * sens
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractVector, sens, ::KA.CPU)
    sig .= vec(transpose(Mxy) * sens)
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractVector, sens, ::KA.GPU)
    sig .= vec(sum(reshape(Mxy, :, 1) .* sens; dims=1))
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractMatrix, sens, backend::KA.GPU)
    for sample in axes(Mxy, 2)
        combine_coil_signal!(
            @view(sig[sample, :]), @view(Mxy[:, sample]), sens, backend,
        )
    end
    return nothing
end

combine_coil_signal!(sig, Mxy, sens) =
    combine_coil_signal!(sig, Mxy, sens, KA.get_backend(sig))

KomaMRIBase.get_sens(receiver, x, y, z, ::KA.Backend) =
    get_sens(receiver, x, y, z)

struct CoilSensitivities{V,I}
    values::V
    interpolators::I
end

prealloc_sensitivities(receiver::UniformCoilSens, _) = receiver
prealloc_sensitivities(receiver::AbstractRFReceiveSystem, obj) =
    CoilSensitivities(get_sens(receiver, obj.x, obj.y, obj.z), nothing)
prealloc_sensitivities(receiver::ArbitraryCoilSens, obj) = CoilSensitivities(
    get_sens(receiver, obj.x, obj.y, obj.z),
    KomaMRIBase.coil_interpolators(receiver),
)

update_sensitivities!(::UniformCoilSens, _, _, _) = nothing
update_sensitivities!(::CoilSensitivities, _, _, ::NoMotion) = nothing
function update_sensitivities!(
    sensitivities::CoilSensitivities, receiver::AbstractRFReceiveSystem,
    (x, y, z), ::Union{Motion,MotionList},
)
    KomaMRIBase.get_sens!(sensitivities.values, receiver, x, y, z)
    return nothing
end
function update_sensitivities!(
    sensitivities::CoilSensitivities, receiver::ArbitraryCoilSens,
    (x, y, z), ::Union{Motion,MotionList},
)
    KomaMRIBase.get_sens!(
        sensitivities.values, receiver, x, y, z, sensitivities.interpolators,
    )
    return nothing
end

function acquire_signal!(sig, Mxy::AbstractVector, ::UniformCoilSens, _)
    sig .= sum(Mxy)
    return nothing
end

@inline acquire_signal!(sig, sample, Mxy, sensitivities, coords) =
    acquire_signal!(@view(sig[sample, :]), Mxy, sensitivities, coords)

@inline _sample_counter(sig) = firstindex(sig, 1)

# Fills row `sample` of `sig` with `acquire!(row)` when `adc` is on. With a host index the
# row is written in place and skipped when off. An index that is only known at run time
# (e.g. traced by Reactant) cannot branch or be viewed, so the row is stored with a mask.
@inline function store_sample!(acquire!, sig, sample::Integer, adc::Bool)
    adc && acquire!(@view(sig[sample, ntuple(_ -> Colon(), ndims(sig) - 1)...]))
    return nothing
end
@inline function store_sample!(acquire!, sig, sample, adc)
    dims = ntuple(_ -> Colon(), ndims(sig) - 1)
    row = similar(sig, size(sig)[2:end])
    acquire!(row)
    sig[sample, dims...] = adc === true ? row : ifelse.(adc, row, sig[sample, dims...])
    return nothing
end

@inline acquire_signal!(sig, sample, adc, M, sys, prealloc, p, seq, i, _sim_method=nothing) =
    store_sample!(sig, sample, adc) do row
        coords = spin_coordinates!(
            prealloc.coordinates, p.motion, p.x, p.y, p.z, seq.t[i + 1],
        )
        update_sensitivities!(prealloc.sens, sys.receiver, coords, p.motion)
        acquire_signal!(row, M.xy, prealloc.sens, coords)
    end

@inline acquire_signal!(sig, sample, adc, M, sys, prealloc, p, seq, i, block_time, coords) =
    store_sample!(sig, sample, adc) do row
        @. prealloc.M.xy = exp(-block_time / p.T2) * M.xy * cis(prealloc.ϕ)
        outflow_spin_reset!(prealloc.M.xy, spin_reset_state(prealloc))
        update_sensitivities!(prealloc.sens, sys.receiver, coords, p.motion)
        acquire_signal!(row, prealloc.M.xy, prealloc.sens, coords)
    end

function acquire_signal!(sig, Mxy::AbstractMatrix, ::UniformCoilSens, _)
    sig .= vec(sum(Mxy; dims=1))
    return nothing
end

function acquire_signal!(sig, Mxy, sensitivities::CoilSensitivities, _)
    combine_coil_signal!(sig, Mxy, sensitivities.values)
    return nothing
end

function acquire_signal!(sig, Mxy, receiver, (x, y, z))
    combine_coil_signal!(
        sig, Mxy, get_sens(receiver, x, y, z, KA.get_backend(sig)),
    )
    return nothing
end
