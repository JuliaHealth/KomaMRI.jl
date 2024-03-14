using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:files in ti.tags) #verbose=true

@testitem "Files" tags=[:files] begin
    using Suppressor

    # Test Pulseq
    @testset "Pulseq" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000
        @test seq.DEF["signature"] == "67ebeffe6afdf0c393834101c14f3990"

        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000
        @test seq.DEF["signature"] == "efc5eb7dbaa82aba627a31ff689c8649"

        seq = @suppress read_seq(path*"/test_files/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001
        @test seq.DEF["signature"] == "f291a24409c3e8de01ddb93e124d9ff2"

        seq = @suppress read_seq(path*"/test_files/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001
        @test seq.DEF["signature"] == "e827cfff4436b65a6341a4fa0f6deb07"

        # Test Pulseq compression-decompression
        shape = ones(100)
        num_samples, compressed_data = KomaMRIFiles.compress_shape(shape)
        shape2 = KomaMRIFiles.decompress_shape(num_samples, compressed_data)
        @test shape == shape2
    end
    # Test JEMRIS
    @testset "JEMRIS" begin
        path = @__DIR__
        obj = read_phantom_jemris(path*"/test_files/column1d.h5")
        @test obj.name == "column1d.h5"
    end
    # Test JEMRIS
    @testset "MRiLab" begin
        path = @__DIR__
        filename = path * "/test_files/brain_mrilab.mat"
        FRange_filename = path * "/test_files/FRange.mat" #Slab within slice thickness
        obj = read_phantom_MRiLab(filename; FRange_filename)
        @test obj.name == "brain_mrilab.mat"
    end
end

@testitem "Pulseq_Read_Comparison" tags=[:files] begin
    using KomaMRIBase, MAT
    TOLERANCE = 1e-6

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

    @testset "read" begin

        files = ["gr-trapezoidal"; "gr-uniformly-shaped"; "gr-time-shaped"; "rf-pulse"; "rf-uniformly-shaped"; "fid"]
        path = @__DIR__
        for filename in files:

            # Koma Read
            seq = read_seq("$(path)/test_files/pulseq_read_comparison/$(filename).seq")
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
            pulseq_read = matread("$(path)/test_files/pulseq_read_comparison/$(filename).mat")
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

            @test (all([all(.≈(koma_read["rfTimes"][i], pulseq_read["rfTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))
            @test (all([all(.≈(abs.(koma_read["rfAmplitudes"][i] .- pulseq_read["rfAmplitudes"][i]), 0, atol=TOLERANCE)) for i in 1:length(pulseq_read["rfBlocks"])]))
            @test (all([all(.≈(koma_read["gxTimes"][i], pulseq_read["gxTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))
            @test (all([all(.≈(koma_read["gxAmplitudes"][i], pulseq_read["gxAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gxBlocks"])]))
            @test (all([all(.≈(koma_read["gyTimes"][i], pulseq_read["gyTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))
            @test (all([all(.≈(koma_read["gyAmplitudes"][i], pulseq_read["gyAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gyBlocks"])]))
            @test (all([all(.≈(koma_read["gzTimes"][i], pulseq_read["gzTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))
            @test (all([all(.≈(koma_read["gzAmplitudes"][i], pulseq_read["gzAmplitudes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["gzBlocks"])]))
            @test (all([all(.≈(koma_read["adcTimes"][i], pulseq_read["adcTimes"][i], atol=TOLERANCE)) for i in 1:length(pulseq_read["adcBlocks"])]))

        end

    end

end
