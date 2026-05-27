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

function _flow_spin_weights(spin_span, prototype::AbstractVector{T}) where {T}
   w = similar(prototype)
   return spin_indicator(spin_span, w)
end

_flow_reset_mask(mask::AbstractVector) = @view mask[lastindex(mask):lastindex(mask)]
_flow_reset_mask(mask::AbstractMatrix) = selectdim(mask, 2, size(mask, 2))

flow_replace_mag(replace_by::AbstractVector) = replace_by
flow_replace_mag(replace_by::AbstractMatrix) = selectdim(replace_by, 2, size(replace_by, 2))
flow_replace_mag(replace_by::Number) = replace_by

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
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_curve)
   prototype = ndims(spin_state_matrix) == 1 ? spin_state_matrix : selectdim(spin_state_matrix, 2, 1)
   w = _flow_spin_weights(spin_span, prototype)
   mask = get_mask(action.spin_reset, ts)
   indicator = reshape(w, :, 1)
   weighted_mask = mask .* indicator
   spin_state_matrix .*= (1 .- weighted_mask)
   spin_state_matrix .+= replace_by .* weighted_mask
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
   ts = KomaMRIBase.unit_time(init_time(t, seq_t, add_t0), time_curve)
   w = _flow_spin_weights(spin_span, M.z)
   mask = _flow_reset_mask(get_mask(action.spin_reset, ts))
   ρ = flow_replace_mag(replace_by)
   M.xy .*= (1 .- mask .* w)
   M.z  .*= (1 .- mask .* w)
   M.xy .+= zero(eltype(M.z)) .* mask .* w
   M.z  .+= ρ .* mask .* w
   return nothing
end

function init_time(t, seq_t::AbstractArray, add_t0)
   t1 = @view(seq_t[1])
   return add_t0 ? [t1 (t1 .+ t)] : t1 .+ t
end
function init_time(t, seq_t, add_t0)
   return t
end

function get_mask(spin_reset, t)
   itp  = KomaMRIBase.interpolate(spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Next}()), Val(size(spin_reset, 1)), t)
   return KomaMRIBase.resample(itp, t)
end