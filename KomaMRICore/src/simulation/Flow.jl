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
   # Obtain mask
   itp  = KomaMRIBase.interpolate(action.spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(action.spin_reset, 1)), t)
   mask = KomaMRIBase.resample(itp, ts)
   mask .= (cumsum(mask; dims=2) .== 0)
   # Modify spin state: reset and replace by initial value
   spin_state_matrix .*= mask
   spin_state_matrix .+= replace_by .* (1 .- mask)
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
   # Obtain mask
   itp  = KomaMRIBase.interpolate(action.spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(action.spin_reset, 1)), t)
   mask = KomaMRIBase.resample(itp, ts)
   mask .= (cumsum(mask; dims=2) .== 0)
   mask = @view(mask[:, end])
   # Modify spin state: reset and replace by initial value
   M.xy .*= mask
   M.z  .*= mask
   M.xy .+= 0          .* (1 .- mask)
   M.z  .+= replace_by .* (1 .- mask)
   return nothing
end

init_time(t, seq_t, add_t0) = t
init_time(t, seq_t::AbstractArray, add_t0) = begin
   t1 = @view(seq_t[1])
   return add_t0 ? [t1 (t1 .+ t)] : t1 .+ t
end
