"""
    seq = build_label(type, label, value; sys=Scanner())

Return a one-block `Sequence` with a Pulseq-style label extension. See
`make_label` for label arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the label block.
"""
function build_label(type, label, value; sys=Scanner())
    seq = Sequence(sys)
    addblock!(seq, make_label(type, label, value))
    return seq
end

"""
    event = make_label(type, label, value)

Return a Pulseq-style label extension.

# Arguments
- `type`: Label operation, `:SET` or `:INC`.
- `label`: ADC label name, for example `:LIN`, `:SLC`, `:REP`, or `:REF`.
- `value`: Integer label value; `Bool` values map to `0` or `1`.

# Returns
- `event`: `LabelSet` or `LabelInc` extension event.
"""
function make_label(type, label, value)
    return label_event(type)(label_value(value), label_name(label))
end

label_event(type) = error("Label type must be a Symbol, for example :SET or :INC.")
label_event(type::Symbol) = label_event(Val(type))
label_event(::Val{:SET}) = KomaMRIBase.LabelSet
label_event(::Val{:INC}) = KomaMRIBase.LabelInc
label_event(::Val{type}) where {type} =
    error("Unsupported label type `:$type`; use :SET or :INC.")

label_name(label) = error("Label name must be a Symbol, for example :LIN or :SLC.")
function label_name(label::Symbol)
    label in KomaMRIBase.ADC_LABEL_NAMES || error("Unsupported ADC label `$label`.")
    return String(label)
end

label_value(value) = error("Label value must be integer-valued.")
label_value(value::Bool) = Int(value)
function label_value(value::Real)
    isinteger(value) || error("Label value must be integer-valued.")
    return Int(value)
end
