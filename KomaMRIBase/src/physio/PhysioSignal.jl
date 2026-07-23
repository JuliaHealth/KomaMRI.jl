"""
    AbstractPhysioSignal

Supertype for physiological signals used to resolve sequence triggers.
"""
abstract type AbstractPhysioSignal end

_is_trigger(::Extension) = false
_is_trigger(trigger::Trigger) = trigger.type == 2

has_trigger(seq::Sequence) = any(_is_trigger, Iterators.flatten(seq.EXT))
