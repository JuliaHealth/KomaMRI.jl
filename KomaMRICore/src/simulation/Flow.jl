"""
    reset_magnetization!
"""
function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, motion, t, ρ) where {T<:Real}
   return nothing
end

function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, motion::MotionList{T}, t, ρ) where {T<:Real}
   for m in motion.motions
      t_unit = KomaMRIBase.unit_time(t, m.time)
      idx = KomaMRIBase.get_indexing_range(m.spins)
      reset_magnetization!(@view(M[idx]), @view(Mxy[idx, :]), m.action, t_unit, @view(ρ[idx]))
   end
   return nothing
end

function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, action::FlowPath{T}, t, ρ) where {T<:Real}
   itp = KomaMRIBase.interpolate(action.spin_reset, KomaMRIBase.Gridded(KomaMRIBase.Constant{KomaMRIBase.Previous}()), Val(size(action.spin_reset, 1)))
   flags = KomaMRIBase.resample(itp, t)
   reset = vec(any(flags .> 0; dims=2))
   flags = (cumsum(flags; dims=2) .== 0)
   Mxy .*= flags
   M.z[reset] = ρ[reset]
   return nothing
end