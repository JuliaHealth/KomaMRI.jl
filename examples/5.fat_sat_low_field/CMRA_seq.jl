using MAT

# Sequence parameters
# TR = 4.9e-3   #RF Normal
TR = 5.3e-3    #RF Low SAR
iNAV_lines = 6
im_flip_angle = 110
#Pre-pulses
Tfatsat = 26.624e-3
FatSat_flip_angle = 180
T2prep_duration = 50e-3
# Other parameters
TE = TR / 2 #bSSFP condition
Tadc = 1e-6
RR = 1.0 #1 [s]
Trf = 500e-6  #1 [ms]
B1 = 1 / (360*γ*Trf)
im_segments = 20
iNAV_flip_angle = 3.2
number_dummy_heart_beats = 3

function FatSat(α, Δf; sample=false)
    RF_wf = matread("./examples/5.fat_sat_low_field/GAUSS5120.mat")["B1"]
    seq = Sequence()
    seq += Grad(-8e-3, 3000e-6, 500e-6) #Spoiler1
    seq += RF(RF_wf, Tfatsat, Δf)
    α_ref = get_flip_angles(seq)[3]
    seq *= (α/α_ref+0im)
    if sample
        seq += ADC(1, 1e-6)
    end
    seq += Grad(8e-3, 3000e-6, 500e-6) #Spoiler2
    if sample
        seq += ADC(1, 1e-6)
    end
    return seq
end

function T2prep(TE; sample=false)
    seq = Sequence()
    seq += RF(90 * B1, Trf)
    seq += sample ? ADC(20, TE/2 - 1.5Trf) : Delay(TE/2 - 1.5Trf)
    seq += RF(180im * B1 / 2, Trf*2)
    seq += sample ? ADC(20, TE/2 - 1.5Trf) : Delay(TE/2 - 1.5Trf)
    seq += RF(-90 * B1, Trf)
    seq += Grad(8e-3, 6000e-6, 600e-6)
    if sample
        seq += ADC(1, 1e-6)
    end
    return seq
end

function bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=false)
    k = 0
    seq = Sequence()
    for i = 0 : iNAV_lines + im_segments - 1
        if iNAV_lines != 0
            m = (im_flip_angle - iNAV_flip_angle) / iNAV_lines
            α = min( m * i + iNAV_flip_angle, im_flip_angle ) * (-1)^k
        else
            α = im_flip_angle * (-1)^k
        end
        seq += RF(α * B1, Trf)
        if i < iNAV_lines && !sample
            seq += Delay(TR - Trf)
        else
            seq += Delay(TE - Trf/2 - Tadc/2)
            seq += ADC(1, Tadc)
            seq += Delay(TR - TE - Tadc/2 - Trf/2)
        end
        k += 1
    end
    return seq
end

function CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, 
        im_segments, iNAV_flip_angle, im_flip_angle,
        T2prep_duration=50e-3, FatSat_flip_angle=180, RR=1.0; sample_recovery=ones(Bool,number_dummy_heart_beats+1))
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        t2p = T2prep(T2prep_duration; sample=sample_recovery[hb+1])
        fatsat = FatSat(FatSat_flip_angle, fat_freq; sample=sample_recovery[hb+1])
        bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=sample_recovery[hb+1])
        RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
        seq += t2p
        seq += fatsat
        seq += bssfp
        seq += sample_recovery[hb+1] ? ADC(80, RRdelay) : Delay(RRdelay) #Sampling recovery curve
    end
    return seq
end