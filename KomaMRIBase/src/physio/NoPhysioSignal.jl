"""
    NoPhysioSignal()

Marker used when no physiological signal should modify sequence timing.
"""
struct NoPhysioSignal <: AbstractPhysioSignal end

resolve_triggers(seq::Sequence, ::NoPhysioSignal) = seq
