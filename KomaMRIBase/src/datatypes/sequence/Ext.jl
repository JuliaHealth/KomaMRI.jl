export Extension, LabelInc, LabelSet

abstract type Extension end

mutable struct LabelInc <: Extension
  labelstring::String
  labelvalue::Int
end

mutable struct LabelSet <: Extension
  labelstring::String
  labelvalue::Int
end

#mutable struct Trigger <: Extension 
#end