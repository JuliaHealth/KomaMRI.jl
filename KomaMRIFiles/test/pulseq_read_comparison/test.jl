using KomaMRI, MAT, PrettyTables

# Change the filename here
filename = "fid"    # fid, gre, epi

# Auxiliar functions
function get_theo_t_aux(rf::RF)
    return KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
end
function get_theo_A_aux(rf::RF)
    return γ*KomaMRIBase.get_theo_A(rf; off_val=NaN, max_rf_samples=Inf)[2:end]
end
function get_theo_t_aux(gr::Grad)
    return KomaMRIBase.get_theo_t(gr)[2:end]
end
function get_theo_A_aux(gr::Grad)
    return γ*KomaMRIBase.get_theo_A(gr; off_val=NaN)[2:end]
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
    "rfTimes"      => [get_theo_t_aux(seq.RF[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
    "rfAmplitudes" => [get_theo_A_aux(seq.RF[1,i]) for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
    "gxBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gxTimes"      => [get_theo_t_aux(seq.GR[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gxAmplitudes" => [get_theo_A_aux(seq.GR[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
    "gyBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gyTimes"      => [get_theo_t_aux(seq.GR[2,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gyAmplitudes" => [get_theo_A_aux(seq.GR[2,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
    "gzBlocks"     => [i for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
    "gzTimes"      => [get_theo_t_aux(seq.GR[3,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
    "gzAmplitudes" => [get_theo_A_aux(seq.GR[3,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
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
println("rfTimes: $(all([all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=1e-12)) for i in 1:length(pulseq_read["rfBlocks"])]))")
println("rfAmplitudes: $(all([all(.≈(koma_read["rfAmplitudes"][i], pulseq_read["rfAmplitudes"][i], atol=1e-12)) for i in 1:length(pulseq_read["rfBlocks"])]))")
println("gxTimes: $(all([all(.≈(koma_read["gxTimes"][i], pulseq_read["gxTimes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gxBlocks"])]))")
println("gxAmplitudes: $(all([all(.≈(koma_read["gxAmplitudes"][i], pulseq_read["gxAmplitudes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gxBlocks"])]))")
println("gyTimes: $(all([all(.≈(koma_read["gyTimes"][i], pulseq_read["gyTimes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gyBlocks"])]))")
println("gyAmplitudes: $(all([all(.≈(koma_read["gyAmplitudes"][i], pulseq_read["gyAmplitudes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gyBlocks"])]))")
println("gzTimes: $(all([all(.≈(koma_read["gzTimes"][i], pulseq_read["gzTimes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gzBlocks"])]))")
println("gzAmplitudes: $(all([all(.≈(koma_read["gzAmplitudes"][i], pulseq_read["gzAmplitudes"][i], atol=1e-12)) for i in 1:length(pulseq_read["gzBlocks"])]))")
println("adcTimes: $(all([all(.≈(koma_read["adcTimes"][i], pulseq_read["adcTimes"][i], atol=1e-12)) for i in 1:length(pulseq_read["adcBlocks"])]))")

# Just for checking things
#for i in 1:length(pulseq_read["rfBlocks"])
#    println(all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=1e-12)))
#    #println(all(koma_read["rfTimes"][i] .== pulseq_read["rfTimes"][i]))
#end
#for i in 1:length(pulseq_read["rfBlocks"])
#    println("rfTimes: Block $(pulseq_read["rfBlocks"][i])")
#    data = [koma_read["rfTimes"][i] pulseq_read["rfTimes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
#for i in 1:length(pulseq_read["rfBlocks"])
#    println("rfAmplitudes: Block $(pulseq_read["rfBlocks"][i])")
#    data = [koma_read["rfAmplitudes"][i] pulseq_read["rfAmplitudes"][i]]
#    pretty_table(data; header=["koma", "pulseq"])
#end
