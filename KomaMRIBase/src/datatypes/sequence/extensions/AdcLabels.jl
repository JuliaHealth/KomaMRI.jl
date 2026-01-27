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
        SMS::Int=0) = new(LIN, PAR, SLC, SEG, REP, AVG, SET, ECO, PHS, NAV, REV, SMS)
end
  
import Base.show
function Base.show(io::IO, label::AdcLabels)
    label_string = "AdcLabels[ "
    for field in fieldnames(AdcLabels)
        label_string = label_string * string(field)*" = "* string(getfield(label, field))*" | "
    end
    print(io, label_string[1:end-2] * "]")
end