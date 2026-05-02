struct AdcLabels
    # Counters
    LIN::Int   # Copied to Siemens MDH
    PAR::Int   # Copied to Siemens MDH
    SLC::Int   # Copied to Siemens MDH
    SEG::Int   # Copied to Siemens MDH
    REP::Int   # Copied to Siemens MDH
    AVG::Int   # Copied to Siemens MDH
    SET::Int   # Copied to Siemens MDH
    ECO::Int   # Copied to Siemens MDH
    PHS::Int   # Copied to Siemens MDH
    ACQ::Int   # Copied to Siemens MDH
    TRID::Int  # TR segment id for GE and some other interpreters

    # Flags
    NAV::Int    # Copied to Siemens MDH
    REV::Int    # Copied to Siemens MDH
    SMS::Int    # Copied to Siemens MDH
    REF::Int    # Parallel imaging
    IMA::Int    # Parallel imaging
    OFF::Int    # Offline data; negates Siemens ONLINE MDH flag
    NOISE::Int  # Noise-adjust scan for parallel imaging

    # Control
    PMC::Int    # MoCo/PMC blocks for prospective motion correction
    NOROT::Int  # Ignore UI FOV rotation
    NOPOS::Int  # Ignore UI FOV position
    NOSCL::Int  # Ignore UI FOV scaling
    ONCE::Int   # 0: every repeat, 1: first repeat, 2: last repeat

    AdcLabels(
        LIN::Int=0,
        PAR::Int=0,
        SLC::Int=0,
        SEG::Int=0,
        REP::Int=0,
        AVG::Int=0,
        SET::Int=0,
        ECO::Int=0,
        PHS::Int=0,
        ACQ::Int=0,
        TRID::Int=0,
        NAV::Int=0,
        REV::Int=0,
        SMS::Int=0,
        REF::Int=0,
        IMA::Int=0,
        OFF::Int=0,
        NOISE::Int=0,
        PMC::Int=0,
        NOROT::Int=0,
        NOPOS::Int=0,
        NOSCL::Int=0,
        ONCE::Int=0) = new(LIN, PAR, SLC, SEG, REP, AVG, SET, ECO, PHS, ACQ, TRID, NAV, REV, SMS, REF, IMA, OFF, NOISE, PMC, NOROT, NOPOS, NOSCL, ONCE)
end

const ADC_LABEL_NAMES = fieldnames(AdcLabels)

function _adc_label_with(label::AdcLabels, field::Symbol, value::Integer)
    idx = findfirst(==(field), ADC_LABEL_NAMES)
    isnothing(idx) && error("Unsupported ADC label `$field`.")
    values = ntuple(i -> i == idx ? Int(value) : getfield(label, i), Val(length(ADC_LABEL_NAMES)))
    return AdcLabels(values...)
end
  
import Base.show
function Base.show(io::IO, label::AdcLabels)
    label_string = "AdcLabels[ "
    for field in fieldnames(AdcLabels)
        label_string = label_string * string(field)*" = "* string(getfield(label, field))*" | "
    end
    print(io, label_string[1:end-2] * "]")
end
