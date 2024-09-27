function outflow_spin_reset!(args...; kwargs...)
   return nothing
end

function outflow_spin_reset!(
   spin_state_matrix, t, motion::MotionList; 
   replace_by=0, seq_t=0, add_t0=false
)
   for m in motion.motions   
      outflow_spin_reset!(
         spin_state_matrix, t, m.action, m.time, m.spins; 
         replace_by=replace_by, seq_t=seq_t, add_t0=add_t0
      )
   end
   return nothing
end

function outflow_spin_reset!(
    spin_state_matrix::AbstractArray,
    t,
    action::FlowPath,
    time_span,
    spin_span;
    replace_by=0,
    seq_t=0,
    add_t0=false,
)
   # Initialize time: add t0 and normalize
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_span)
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
    time_span,
    spin_span;
    replace_by=0,
    seq_t=0,
    add_t0=false,
)
   # Initialize time: add t0 and normalize
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_span)
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

function get_mask(spin_reset, t::Real)
   itp  = KomaMRIBase.interpolate(spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(spin_reset, 1)), t)
   return KomaMRIBase.resample(itp, t)
end
function get_mask(spin_reset, t::AbstractArray)
   itp  = KomaMRIBase.interpolate(spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(spin_reset, 1)), t)
   mask = KomaMRIBase.resample(itp, t)
   mask .= (cumsum(mask; dims=2) .== 1)
   return mask
end