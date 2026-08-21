spin_coordinates(motion, x, y, z, t) = get_spin_coords(motion, x, y, z, t)
spin_coordinates(::NoMotion, x, y, z, t) = x, y, z

struct CycleRemapWorkspace{I,M}
    source::I
    scratch::M
end

cycle_remap_indices(::Colon, N) = Base.OneTo(N)
cycle_remap_indices(indices, N) = indices

function cycle_remap_sources(obj::Phantom, motion::Motion)
    affected = collect(cycle_remap_indices(KomaMRIBase.get_indexing_range(motion.spins), length(obj)))
    map = motion.action.cycle_map
    length(affected) == length(map) ||
        throw(DimensionMismatch("cycle_map length must equal the number of affected phantom particles."))
    source = collect(eachindex(obj.ρ))
    source[affected] .= affected[map]
    return source
end

function cycle_remap_break_indices(seqd, motion)
    breaks = Int[]
    for t in KomaMRIBase.cycle_remap_times(motion, last(seqd.t))
        t < last(seqd.t) || continue
        i = searchsortedlast(seqd.t, t)
        i >= firstindex(seqd.t) && seqd.t[i] == t ||
            error("Cycle boundary $t is missing from the simulation sampling grid.")
        push!(breaks, i)
    end
    return breaks
end

@kernel function cycle_remap_kernel!(M_xy_new, M_z_new, @Const(M_xy), @Const(M_z), @Const(source))
    i = @index(Global)
    @inbounds begin
        M_xy_new[i] = M_xy[source[i]]
        M_z_new[i] = M_z[source[i]]
    end
end

function cycle_remap_workspace(source, M::Mag, backend::KA.GPU)
    return CycleRemapWorkspace(
        gpu(source, backend),
        Mag(similar(M.xy), similar(M.z)),
    )
end

function remap_magnetization!(M::Mag, workspace::CycleRemapWorkspace, backend::KA.GPU, groupsize)
    cycle_remap_kernel!(backend, groupsize)(
        workspace.scratch.xy,
        workspace.scratch.z,
        M.xy,
        M.z,
        workspace.source;
        ndrange=length(workspace.source),
    )
    KA.synchronize(backend)
    M.xy, workspace.scratch.xy = workspace.scratch.xy, M.xy
    M.z, workspace.scratch.z = workspace.scratch.z, M.z
    return nothing
end

outflow_spin_reset_at!(spin_state, t, i, motion; replace_by=0) =
   outflow_spin_reset!(spin_state, t[i, :], motion; replace_by)
outflow_spin_reset_at!(spin_state, t, i, ::NoMotion; replace_by=0) = nothing

function outflow_spin_reset!(args...; kwargs...)
   return nothing
end

function outflow_spin_reset!(spin_state_matrix, t, ml::MotionList; replace_by=0, seq_t=0, add_t0=false)
   for m in ml.motions   
      outflow_spin_reset!(spin_state_matrix, t, m; replace_by=replace_by, seq_t=seq_t, add_t0=add_t0)
   end
   return nothing
end

function outflow_spin_reset!(spin_state_matrix, t, m::Motion; replace_by=0, seq_t=0, add_t0=false) 
   outflow_spin_reset!(spin_state_matrix, t, m.action, m.time, m.spins; replace_by=replace_by, seq_t=seq_t, add_t0=add_t0)
   return nothing
end

function outflow_spin_reset!(
    spin_state_matrix::AbstractArray,
    t,
    action::FlowPath,
    time_curve,
    spin_span;
    replace_by=0,
    seq_t=0,
    add_t0=false,
)
   # Initialize time: add t0 and normalize
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_curve)
   # Get spin state range affected by the spin span
   idx = KomaMRIBase.get_indexing_range(spin_span)
   spin_state_matrix = @view(spin_state_matrix[idx, :])
   replace_by = replace_view(replace_by, idx)
   # Obtain mask
   mask = get_mask(action.spin_reset, ts)
   # Modify spin state: reset and replace by initial value
   spin_state_matrix .*= (1 .- mask)
   spin_state_matrix .+= replace_by .* mask
   return nothing
end

function outflow_spin_reset!(
    M::Mag,
    t,
    action::FlowPath,
    time_curve,
    spin_span;
    replace_by=0,
    seq_t=0,
    add_t0=false,
)
   # Initialize time: add t0 and normalize
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_curve)
   # Get spin state range affected by the spin span
   idx = KomaMRIBase.get_indexing_range(spin_span)
   M = @view(M[idx])
   replace_by = replace_view(replace_by, idx)
   # Obtain mask
   mask = get_mask(action.spin_reset, ts)
   mask = @view(mask[:, end])
   # Modify spin state: reset and replace by initial value
   M.xy .*= (1 .- mask)
   M.z  .*= (1 .- mask)
   M.xy .+= 0          .* mask
   M.z  .+= replace_by .* mask
   return nothing
end

function init_time(t, seq_t::AbstractArray, add_t0)
   t1 = @view(seq_t[1])
   return add_t0 ? [t1 (t1 .+ t)] : t1 .+ t
end
function init_time(t, seq_t, add_t0)
   return t
end

function replace_view(replace_by::AbstractArray, idx)
   return @view(replace_by[idx])
end
function replace_view(replace_by, idx)
   return replace_by
end

function get_mask(spin_reset, t)
   itp  = KomaMRIBase.interpolate(spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Next}()), Val(size(spin_reset, 1)), t)
   return KomaMRIBase.resample(itp, t)
end
