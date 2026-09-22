spin_coordinates(motion, x, y, z, t) = get_spin_coords(motion, x, y, z, t)
spin_coordinates(::NoMotion, x, y, z, t) = x, y, z

struct NoSpinReset end

struct SpinResetState{I,B,C,S}
   # Earliest simulation step with a reset flag, and its running block state.
   first_reset::I
   has_reset::B
   # Source motion-column labels and their sampled value at each simulation step.
   reset_column_ids::C
   sampled_reset_columns::S
end

Base.view(::NoSpinReset, _) = NoSpinReset()
function Base.view(state::SpinResetState, idx)
   return SpinResetState(
      @view(state.first_reset[idx]),
      @view(state.has_reset[idx]),
      state.reset_column_ids,
      state.sampled_reset_columns,
   )
end

prealloc_view(state::SpinResetState, idx) = view(state, idx)
spin_reset_state(::PreallocResult) = NoSpinReset()
spin_reset_state(prealloc::BlochSimplePrealloc) = prealloc.spin_reset
spin_reset_state(prealloc::BlochCPUPrealloc) = prealloc.spin_reset
spin_reset_state(prealloc::BlochGPUPrealloc) = prealloc.spin_reset
spin_reset_state(prealloc::BlochMagnusCPUPrealloc) = prealloc.spin_reset

first_reset_steps(state::SpinResetState) = state.first_reset
first_reset_steps(::NoSpinReset) = nothing

reset_sample_count(_) = 0
reset_sample_count(action::FlowPath) = size(action.spin_reset, 2)
reset_sample_count(motion::Motion) = reset_sample_count(motion.action)
reset_sample_count(motion::MotionList) =
   maximum(reset_sample_count(m) for m in motion.motions; init=0)

function prealloc_spin_reset(motion, backend, nspins, max_block_length)
   max_columns = reset_sample_count(motion)
   iszero(max_columns) && return NoSpinReset()

   reset_column_ids = KA.zeros(backend, UInt32, max_columns)
   copyto!(reset_column_ids, collect(UInt32, 1:max_columns))
   return SpinResetState(
      KA.zeros(backend, UInt32, nspins),
      KA.zeros(backend, Bool, nspins),
      reset_column_ids,
      KA.zeros(backend, UInt32, max_block_length),
   )
end

KA.@kernel unsafe_indices=true inbounds=true function first_reset_step_kernel!(
   first_reset,
   @Const(spin_reset),
   @Const(sampled_reset_columns),
   nsteps,
)
   spin = @index(Global, Linear)
   step = UInt32(1)
   first = first_reset[spin]
   while step <= nsteps && step < first
      if spin_reset[spin, sampled_reset_columns[step]]
         first = step
         break
      end
      step += UInt32(1)
   end
   first_reset[spin] = first
end

# Prepares the SpinResetState for a new simulation block
function prepare_spin_reset!(state::SpinResetState, motion, t, backend)
   fill!(state.first_reset, typemax(UInt32))
   fill!(state.has_reset, false)
   accumulate_spin_resets!(state, motion, vec(t), backend)
   return nothing
end
prepare_spin_reset!(::NoSpinReset, motion, t, backend) = nothing

# Accumulates the earliest simulation step with a reset flag for each spin
function accumulate_spin_resets!(state, motion::MotionList, t, backend)
   for m in motion.motions
      accumulate_spin_resets!(state, m, t, backend)
   end
   return nothing
end
function accumulate_spin_resets!(state, motion::Motion, t, backend)
   accumulate_spin_resets!(state, motion.action, motion.time, motion.spins, t, backend)
   return nothing
end
accumulate_spin_resets!(state, action, time, spins, t, backend) = nothing
function accumulate_spin_resets!(
   state,
   action::FlowPath,
   time,
   spins,
   t,
   backend,
)
   # Use the existing Constant{Next} convention to map each block time to a
   # column of the FlowPath reset matrix without sampling every spin.
   normalized_time = KomaMRIBase.unit_time(t, time)
   ncolumns = size(action.spin_reset, 2)
   nsteps = length(t)
   reset_column_ids =
      reshape(@view(state.reset_column_ids[1:ncolumns]), 1, ncolumns)
   itp = KomaMRIBase.interpolate(
      reset_column_ids,
      KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Next}()),
      Val(1),
      normalized_time,
   )
   sampled_reset_columns = @view(state.sampled_reset_columns[1:nsteps])
   sampled_reset_columns .= itp.(normalized_time)

   # Retain only the first flagged simulation step for each affected spin.
   idx = KomaMRIBase.get_indexing_range(spins)
   first_reset = @view(state.first_reset[idx])
   first_reset_step_kernel!(backend)(
      first_reset,
      action.spin_reset,
      sampled_reset_columns,
      UInt32(nsteps);
      ndrange=length(first_reset),
   )
   KA.synchronize(backend)
   return nothing
end

# Advances SpinResetState to the current simulation step, setting the has_reset flag
function advance_spin_reset!(state::SpinResetState, step)
   state.has_reset .|= UInt32(step) .>= state.first_reset
   return nothing
end
advance_spin_reset!(::NoSpinReset, step) = nothing

function outflow_spin_reset!(spin_state, state::SpinResetState; replace_by=0)
   @. spin_state = ifelse(state.has_reset, replace_by, spin_state)
   return nothing
end
outflow_spin_reset!(spin_state, ::NoSpinReset; replace_by=0) = nothing

function outflow_spin_reset!(M::Mag, state::SpinResetState; replace_by=0)
   @. M.xy = ifelse(state.has_reset, zero(eltype(M.xy)), M.xy)
   @. M.z = ifelse(state.has_reset, replace_by, M.z)
   return nothing
end

function outflow_spin_reset_at!(
   spin_state, state::SpinResetState, step, t, motion; replace_by=0,
)
   advance_spin_reset!(state, step)
   outflow_spin_reset!(spin_state, state; replace_by)
   return nothing
end
function outflow_spin_reset_at!(
   spin_state, ::NoSpinReset, step, t, motion; replace_by=0,
)
   return nothing
end

@inline spin_first_reset(first_reset, spin) = @inbounds first_reset[spin]
@inline spin_first_reset(::Nothing, spin) = typemax(UInt32)
@inline spin_has_reset(first_reset, step) = UInt32(step) >= first_reset

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
