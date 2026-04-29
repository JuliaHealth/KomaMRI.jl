function round_trip_sequences()
    sequences = Sequence[]

    amp(t) = 10e-6 * sqrt.(π*t) .* sin.(π*t)
    
    # 1. FID sequence
    seq = Sequence()
    seq += RF(1e-6, 100e-6)
    seq += ADC(100, 100e-3, 50e-3)
    seq.DEF["Name"] = "fid"
    push!(sequences, seq)

    # 2. EPI sequence
    seq = PulseDesigner.EPI_example()
    seq.DEF["Name"] = "epi"
    push!(sequences, seq)

    # 3. Uniformly-shaped gradient
    t = 0:0.1:1
    A = amp(t)
    T = 10e-3
    delay, rise, fall = 1e-3, 5e-6 , 5e-6
    first, last = 1e-6, -5e-6
    gr = Grad(A, T, rise, fall, delay, first, last)
    seq = Sequence([gr])
    seq.DEF["Name"] = "gr-uniformly-shaped"
    push!(sequences, seq)

    # 4. Time-shaped gradient
    t = 0:0.1:1
    A = amp(t)
    T = 1e-3 * rand(length(A)-1)
    delay, rise, fall = 1e-3, 5e-6 , 5e-6
    first, last = 1e-6, -5e-6
    gr = Grad(A, T, rise, fall, delay, first, last)
    seq = Sequence([gr])
    seq.DEF["Name"] = "gr-time-shaped"
    push!(sequences, seq)

    # 5. Uniformly-shaped RF with phase offset
    t = 0:0.1:1
    A = 0.5 * amp(t)
    T = 10e-3
    delay = 1e-3
    rf = RF(A, T, 0, delay)
    seq = Sequence(); seq += rf; seq += complex(-1.0) * rf
    seq.DEF["Name"] = "rf-uniformly-shaped"
    push!(sequences, seq)

    # 6. Time-shaped RF with phase offset
    t = 0:0.1:1
    A = 0.5 * amp(t)
    T = 1e-3 * rand(length(A)-1)
    delay = 0.1e-3;
    rf = RF(A, T, 0, delay)
    seq = Sequence(); seq += rf
    seq.DEF["Name"] = "rf-time-shaped"
    push!(sequences, seq)

    # 7. Frequency-modulated RF (adiabatic)
    ξ = 2π
    Bx0 = 1e-5
    Bz0 = 1e-5
    T = π / (2 * ξ) # 90° excitation pulse
    delay = 0.1e-3;
    # Modulation functions: (Ugurbil et al., 1987)
    AM(t) = Bx0 * sin(ξ*t)
    PM(t) = -(γ*Bz0 / ξ) * sin(ξ*t)
    # Implement modulation:
    t = 0:0.01:T
    A = AM.(t) .* exp.(-1im * PM.(t))
    rf = RF(A, T, 0, delay)
    seq = Sequence(); seq += rf
    seq.DEF["Name"] = "rf-frequency-modulated"
    push!(sequences, seq)

    # 8. Combination of Events (see https://juliahealth.org/KomaMRI.jl/dev/explanation/5-seq-events#Combination-of-Events)
    rf = RF(1e-6*[0; -0.1; 0.2; -0.5; 1; -0.5; 0.2; -0.1; 0], 0.5e-3)
    gx = Grad(50*10e-6, 5e-3, 1e-3)
    gy = Grad(50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1], 5e-3, 2e-3)
    gz = Grad(50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1], 1e-3, 1e-3)
    adc = ADC(16, 5e-3, 0)
    seq = Sequence([gx; gy; gz;;], [rf;;], [adc])
    seq.DEF["Name"] = "combination-of-events"
    push!(sequences, seq)

    # 9. Extensions
    seq = Sequence()
    seq += RF(1e-6, 100e-6)
    seq += ADC(100, 100e-3, 50e-3)
    seq.EXT = [[LabelInc(1, "LIN"), LabelSet(1, "ECO")], [LabelSet(1, "ECO")]]
    seq.DEF["Name"] = "extensions"

    # 10. Official Pulseq MATLAB GRE example with RF spoiling
    # https://pulseq.github.io/writeGradientEcho.html
    push!(sequences, seq)
    seq = read_seq(joinpath(@__DIR__, "test_files/pulseq/read_comparison/v1.5/gre.seq"))
    seq.DEF["Name"] = "gre_rfspoiled"
    push!(sequences, seq)

    return sequences
end
