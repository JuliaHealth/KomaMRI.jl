function round_trip_sequences()
    sequences = Sequence[]

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
    A = 10*10e-6 * sqrt.(π*t) .* sin.(π*t)
    T = 10e-3
    delay, rise, fall = 1e-3, 5e-6 , 5e-6
    first, last = 1e-6, -5e-6
    gr = Grad(A, T, rise, fall, delay, first, last)
    seq = Sequence([gr])
    seq.DEF["Name"] = "gr-uniformly-shaped"
    push!(sequences, seq)

    # 4. Time-shaped gradient
    t = 0:0.1:1
    A = 10*10e-6 * sqrt.(π*t) .* sin.(π*t)
    T = 1e-3 * rand(length(A)-1)
    delay, rise, fall = 1e-3, 5e-6 , 5e-6
    first, last = 1e-6, -5e-6
    gr = Grad(A, T, rise, fall, delay, first, last)
    seq = Sequence([gr])
    seq.DEF["Name"] = "gr-time-shaped"
    push!(sequences, seq)

    # 5. Combination of Events (see https://juliahealth.org/KomaMRI.jl/dev/explanation/5-seq-events#Combination-of-Events)
    rf = RF(1e-6*[0; -0.1; 0.2; -0.5; 1; -0.5; 0.2; -0.1; 0], 0.5e-3)
    gx = Grad(50*10e-6, 5e-3, 1e-3)
    gy = Grad(50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1], 5e-3, 2e-3)
    gz = Grad(50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1], 1e-3, 1e-3)
    adc = ADC(16, 5e-3, 0.5e-3)
    seq = Sequence([gx; gy; gz;;], [rf;;], [adc])
    seq.DEF["Name"] = "combination-of-events"
    push!(sequences, seq)

    return sequences
end