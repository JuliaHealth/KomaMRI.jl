spin_coordinates(motion, x, y, z, t) = get_spin_coords(motion, x, y, z, t)
spin_coordinates(::NoMotion, x, y, z, t) = x, y, z

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
   spin_state_matrix .= ifelse.(mask, replace_by, spin_state_matrix)
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
   ts = KomaMRIBase.unit_time(last_time(init_time(t, seq_t, add_t0)), time_curve)
   # Get spin state range affected by the spin span
   idx = KomaMRIBase.get_indexing_range(spin_span)
   M = @view(M[idx])
   replace_by = replace_view(replace_by, idx)
   # Obtain mask
   mask = vec(get_mask(action.spin_reset, ts))
   # Modify spin state: reset and replace by initial value
   M.xy .= ifelse.(mask, zero(eltype(M.xy)), M.xy)
   M.z  .= ifelse.(mask, replace_by, M.z)
   return nothing
end

last_time(t::Real) = t
last_time(t::AbstractArray) = @view t[end:end]

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
   return KomaMRIBase.next_columns(spin_reset, t)
end
