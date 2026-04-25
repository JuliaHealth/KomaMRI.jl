mutable struct AdcLabels
    LIN::Int
    PAR::Int
    SLC::Int
    SEG::Int
    REP::Int
    AVG::Int
    SET::Int
    ECO::Int
    PHS::Int
    NAV::Int
    REV::Int
    SMS::Int
    ACQ::Int
    OFF::Int
    NOISE::Int
    REF::Int
    IMA::Int
    PMC::Int
    NOPOS::Int
    NOROT::Int
    NOSCL::Int
    ONCE::Int
    TRID::Int

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
        NAV::Int=0,
        REV::Int=0,
        SMS::Int=0,
        ACQ::Int=0,
        OFF::Int=0,
        NOISE::Int=0,
        REF::Int=0,
        IMA::Int=0,
        PMC::Int=0,
        NOPOS::Int=0,
        NOROT::Int=0,
        NOSCL::Int=0,
        ONCE::Int=0,
        TRID::Int=0) = new(LIN, PAR, SLC, SEG, REP, AVG, SET, ECO, PHS, NAV, REV, SMS, ACQ, OFF, NOISE, REF, IMA, PMC, NOPOS, NOROT, NOSCL, ONCE, TRID)
end
  
import Base.show
function Base.show(io::IO, label::AdcLabels)
    label_string = "AdcLabels[ "
    for field in fieldnames(AdcLabels)
        label_string = label_string * string(field)*" = "* string(getfield(label, field))*" | "
    end
    print(io, label_string[1:end-2] * "]")
end
