"""
    reset_magnetization!
"""
function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, motion::NoMotion{T}, t::AbstractArray{T}) where {T<:Real}
   return nothing
end

function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, motion::MotionList{T}, t::AbstractArray{T}) where {T<:Real}
   for m in motion.motions
      idx = KomaMRIBase.get_idx(m.spins)
      reset_magnetization!(@view(M[idx]), @view(Mxy[idx, :]), m.action, t)
   end
   return nothing
end

function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, action::KomaMRIBase.AbstractActionSpan{T}, t::AbstractArray{T}) where {T<:Real}
   return nothing
end

function reset_magnetization!(M::Mag{T}, Mxy::AbstractArray{Complex{T}}, action::FlowPath{T}, t::AbstractArray{T}) where {T<:Real}
    itp = interpolate(action.spin_reset, Gridded(Constant{Previous}), Val(size(action.spin_reset, 1)))
    flags = resample(itp, unit_time(t, action.time))
    reset = any(flags; dims=2)
    flags = .!(cumsum(flags; dims=2) .>= 1)
    Mxy .*= flags
    M.z[reset] = p.ρ[reset]
   return nothing
end