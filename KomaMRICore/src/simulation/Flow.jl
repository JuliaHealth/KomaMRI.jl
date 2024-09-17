function outflow_spin_reset!(args...; kwargs...)
   return nothing
end

function outflow_spin_reset!(M, t, motion::MotionList; M0=0, seq_t=0, add_t0=false)
   for m in motion.motions   
      outflow_spin_reset!(M, t, m.action, m.time, m.spins; M0=M0, seq_t=seq_t, add_t0=add_t0)
   end
   return nothing
end

function outflow_spin_reset!(M, t, action::FlowPath, time, spins; M0=0, seq_t=0, add_t0=false)
   t_unit = KomaMRIBase.unit_time(t, time)
   idx = KomaMRIBase.get_indexing_range(spins)
   M = @view(M[idx, :])
   M0 = init_magnetization(M, M0)
   t = init_time(t_unit, seq_t, add_t0)
   itp = KomaMRIBase.interpolate(action.spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(action.spin_reset, 1)))
   mask = KomaMRIBase.resample(itp, t)
   mask .= (cumsum(mask; dims=2) .== 0)
   mask_end = 1 .- vec(any(mask .== 0; dims=2))
   if size(M, 2) > 1
      M .*= mask
      M .+= M0 .* (1 .- mask)
   else
      M .*= mask_end
      M .+= M0 .* (1 .- mask_end)
   end
   return nothing
end

init_time(t, seq_t, add_t0) = t
init_time(t, seq_t::AbstractArray, add_t0) = begin
   t1 = @view(seq_t[1])
   return add_t0 ? [t1 (t1 .+ t)] : t1 .+ t
end

init_magnetization(M, M0) = M0
init_magnetization(M, M0::Real) = begin
   x = similar(M, size(M,1))
   x .= M0
   return x
end

