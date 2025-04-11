export Extension, LabelInc, LabelSet

abstract type Extension end

mutable struct LabelInc <: Extension
  labelvalue::Int
  labelstring::String
end

mutable struct LabelSet <: Extension
  labelvalue::Int
  labelstring::String
end

mutable struct Trigger <: Extension 
  type::Int
  channel::Int
  d1::Float64
  d2::Float64
end