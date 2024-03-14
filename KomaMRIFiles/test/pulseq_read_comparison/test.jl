using KomaMRI, MAT, PrettyTables

TOLERANCE = 1e-6

# Change the filename here
# gr-trapezoidal, gr-uniformly-shaped, gr-time-shaped
# rf-pulse, rf-uniformly-shaped
# fid, spiral
# MSE_test_KomaMRI
filename = "MSE_test_KomaMRI"

# Auxiliar functions
function get_theo_t_aux(rf::RF)
    Namp, Ntim = length(rf.A), length(rf.T)
    if Namp == 1 && Ntim == 1
        return KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
    elseif Namp > 1 && Ntim == 1
        amps = KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
        return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
    end
    return []#KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
end
function get_theo_A_aux(rf::RF)
    Namp, Ntim = length(rf.A), length(rf.T)
    if Namp == 1 && Ntim == 1
        return γ*KomaMRIBase.get_theo_A(rf; off_val=NaN, max_rf_samples=Inf)[2:end]
    elseif Namp > 1 && Ntim == 1
        amps = γ*KomaMRIBase.get_theo_A(rf; off_val=NaN, max_rf_samples=Inf)[2:end]
        return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
    end
    return []#γ*KomaMRIBase.get_theo_A(rf; off_val=NaN, max_rf_samples=Inf)[2:end]
end
function get_theo_t_aux(gr::Grad)
    Namp, Ntim = length(gr.A), length(gr.T)
    if Namp == 1 && Ntim == 1
        return KomaMRIBase.get_theo_t(gr)
    elseif Namp > 1 && Ntim == 1
        return KomaMRIBase.get_theo_t(gr)[2:end]
    elseif Namp > 1 && Ntim > 1
        return KomaMRIBase.get_theo_t(gr)[2:end]
    end
    return []
end
function get_theo_A_aux(gr::Grad)
    Namp, Ntim = length(gr.A), length(gr.T)
    if Namp == 1 && Ntim == 1
        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)
    elseif  Namp > 1 && Ntim == 1
        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)[2:end]
    elseif Namp > 1 && Ntim > 1
        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)[2:end]
    end
    return []
end
function get_theo_t_aux(adc::ADC)
    return (adc.N == 1) ? [adc.T/2] .+ adc.delay : [range(0, adc.T; length=adc.N)...] .+ adc.delay
end

# Koma Read
seq = read_seq("matlab_files/$(filename).seq")
N = length(seq)
T0 = KomaMRIBase.get_block_start_times(seq)
koma_read = Dict(
    "rfBlocks"     => [i for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
    "rfTimes"      => [get_theo_t_aux(seq[i].RF[1]) .+ T0[i] for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
    "rfAmplitudes" => [get_theo_A_aux(seq[i].RF[1]) for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
    "gxBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gxTimes"      => [get_theo_t_aux(seq.GR[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gxAmplitudes" => [get_theo_A_aux(seq.GR[1,i]) for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gyBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gyTimes"      => [get_theo_t_aux(seq.GR[2,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gyAmplitudes" => [get_theo_A_aux(seq.GR[2,i]) for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gzBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
    "gzTimes"      => [get_theo_t_aux(seq.GR[3,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
    "gzAmplitudes" => [get_theo_A_aux(seq.GR[3,i]) for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
    "adcBlocks"    => [i for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
    "adcTimes"     => [get_theo_t_aux(seq.ADC[i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
    )

# Pulseq Read
pulseq_read = matread("matlab_files/$(filename).mat")
pulseq_read["rfBlocks"] = vec(pulseq_read["rfBlocks"])
pulseq_read["rfTimes"]  = vec([vec(blockRfTimes) for blockRfTimes in pulseq_read["rfTimes"]])
pulseq_read["rfAmplitudes"]  = vec([vec(blockRfAmplitudes) for blockRfAmplitudes in pulseq_read["rfAmplitudes"]])
pulseq_read["gxBlocks"] = vec(pulseq_read["gxBlocks"])
pulseq_read["gxTimes"]  = vec([vec(blockGxTimes) for blockGxTimes in pulseq_read["gxTimes"]])
pulseq_read["gxAmplitudes"]  = vec([vec(blockGxAmplitudes) for blockGxAmplitudes in pulseq_read["gxAmplitudes"]])
pulseq_read["gyBlocks"] = vec(pulseq_read["gyBlocks"])
pulseq_read["gyTimes"]  = vec([vec(blockGyTimes) for blockGyTimes in pulseq_read["gyTimes"]])
pulseq_read["gyAmplitudes"]  = vec([vec(blockGyAmplitudes) for blockGyAmplitudes in pulseq_read["gyAmplitudes"]])
pulseq_read["gzBlocks"] = vec(pulseq_read["gzBlocks"])
pulseq_read["gzTimes"]  = vec([vec(blockGzTimes) for blockGzTimes in pulseq_read["gzTimes"]])
pulseq_read["gzAmplitudes"]  = vec([vec(blockGzAmplitudes) for blockGzAmplitudes in pulseq_read["gzAmplitudes"]])
pulseq_read["adcBlocks"] = vec(pulseq_read["adcBlocks"])
pulseq_read["adcTimes"]  = vec([vec(blockAdcTimes) for blockAdcTimes in pulseq_read["adcTimes"]])

# Comparison (true means that all points are equal)
println("rfTimes: $(all([all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))")
println("rfAmplitudes: $(all([all(.≈(abs.(koma_read["rfAmplitudes"][i] .- pulseq_read["rfAmplitudes"][i]), 0, atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))")
println("gxTimes: $(all([all(.≈(koma_read["gxTimes"][i], pulseq_read["gxTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))")
println("gxAmplitudes: $(all([all(.≈(koma_read["gxAmplitudes"][i], pulseq_read["gxAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))")
println("gyTimes: $(all([all(.≈(koma_read["gyTimes"][i], pulseq_read["gyTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))")
println("gyAmplitudes: $(all([all(.≈(koma_read["gyAmplitudes"][i], pulseq_read["gyAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))")
println("gzTimes: $(all([all(.≈(koma_read["gzTimes"][i], pulseq_read["gzTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))")
println("gzAmplitudes: $(all([all(.≈(koma_read["gzAmplitudes"][i], pulseq_read["gzAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))")
println("adcTimes: $(all([all(.≈(koma_read["adcTimes"][i], pulseq_read["adcTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["adcBlocks"])]))")

#for i in 1:length(pulseq_read["rfBlocks"])
#    println(all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=TOLERANCE)))
#    #println(all(koma_read["rfTimes"][i] .== pulseq_read["rfTimes"][i]))
#end

## Just for checking things visually
#
## For RF times
#for i in 1:length(pulseq_read["rfBlocks"])
#    println("rfTimes: Block $(pulseq_read["rfBlocks"][i]), Index $(i)")
#    data = [koma_read["rfTimes"][i] pulseq_read["rfTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("rfTimes: $(all([all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))")

# For RF amplitudes
for i in 1:length(pulseq_read["rfBlocks"])
    println("rfAmplitudes: Block $(pulseq_read["rfBlocks"][i]), Index $(i)")
    data = [koma_read["rfAmplitudes"][i] pulseq_read["rfAmplitudes"][i]]
    pretty_table(data; header=["koma", "pulseq"])
end
println("rfAmplitudes: $(all([all(.≈(koma_read["rfAmplitudes"][i], pulseq_read["rfAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))")

## For GX times
#for i in 1:length(pulseq_read["gxBlocks"])
#    println("gxTimes: Block $(pulseq_read["gxBlocks"][i]), Index $(i)")
#    data = [koma_read["gxTimes"][i] pulseq_read["gxTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gxTimes: $(all([all(.≈(koma_read["gxTimes"][i], pulseq_read["gxTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))")
#
## For GX amplitudes
#for i in 1:length(pulseq_read["gxBlocks"])
#    println("gxAmplitudes: Block $(pulseq_read["gxBlocks"][i]), Index $(i)")
#    data = [koma_read["gxAmplitudes"][i] pulseq_read["gxAmplitudes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gxAmplitudes: $(all([all(.≈(koma_read["gxAmplitudes"][i], pulseq_read["gxAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))")
#
## For GY times
#for i in 1:length(pulseq_read["gyBlocks"])
#    println("gyTimes: Block $(pulseq_read["gyBlocks"][i]), Index $(i)")
#    data = [koma_read["gyTimes"][i] pulseq_read["gyTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gyTimes: $(all([all(.≈(koma_read["gyTimes"][i], pulseq_read["gyTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))")
#
## For GY amplitudes
#for i in 1:length(pulseq_read["gzBlocks"])
#    println("gyAmplitudes: Block $(pulseq_read["gyBlocks"][i]), Index $(i)")
#    data = [koma_read["gyAmplitudes"][i] pulseq_read["gyAmplitudes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gyAmplitudes: $(all([all(.≈(koma_read["gyAmplitudes"][i], pulseq_read["gyAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))")
#
## For GZ times
#for i in 1:length(pulseq_read["gzBlocks"])
#    println("gzTimes: Block $(pulseq_read["gzBlocks"][i]), Index $(i)")
#    data = [koma_read["gzTimes"][i] pulseq_read["gzTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gzTimes: $(all([all(.≈(koma_read["gzTimes"][i], pulseq_read["gzTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))")
#
## For GZ amplitudes
#for i in 1:length(pulseq_read["gzBlocks"])
#    println("gzAmplitudes: Block $(pulseq_read["gzBlocks"][i]), Index $(i)")
#    data = [koma_read["gzAmplitudes"][i] pulseq_read["gzAmplitudes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("gzAmplitudes: $(all([all(.≈(koma_read["gzAmplitudes"][i], pulseq_read["gzAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))")
#
## For ADC
#for i in 1:length(pulseq_read["adcBlocks"])
#    println("adcTimes: Block $(pulseq_read["adcBlocks"][i]), Index $(i)")
#    data = [koma_read["adcTimes"][i] pulseq_read["adcTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#println("adcTimes: $(all([all(.≈(koma_read["adcTimes"][i], pulseq_read["adcTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["adcBlocks"])]))")
#
