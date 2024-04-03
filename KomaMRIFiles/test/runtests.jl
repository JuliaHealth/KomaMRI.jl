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
    using KomaMRIBase, MAT, PrettyTables
    TOLERANCE = 1e-6

    get_theo_A_helper(r::RF) = sum(abs.(r.A)) == 0 ? [0; 0; zeros(2*length(r.A)); 0] : [0; 0; transpose([r.A r.A])[:]; 0]

    # Auxiliar functions
    function get_theo_t_aux(rf::RF)
        Namp, Ntim = length(rf.A), length(rf.T)
        if Namp == 1 && Ntim == 1
            return KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
        elseif Namp > 1 && Ntim == 1
            amps = KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
            return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
        end
        return []
    end
    function get_theo_A_aux(rf::RF)
        Namp, Ntim = length(rf.A), length(rf.T)
        if Namp == 1 && Ntim == 1
            return γ*get_theo_A_helper(rf)[2:end]
        elseif Namp > 1 && Ntim == 1
            amps = γ*get_theo_A_helper(rf)[2:end]
            return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
        end
        return []
    end
    function get_theo_t_aux(gr::Grad)
        if !(gr.A isa Vector) && !(gr.T isa Vector)
            return KomaMRIBase.get_theo_t(gr)
        else
            return KomaMRIBase.get_theo_t(gr)[2:end]
        end
    end
    function get_theo_A_aux(gr::Grad)
        if !(gr.A isa Vector) && !(gr.T isa Vector)
            return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)
        else
            return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)[2:end]
        end
    end
    function get_theo_t_aux(adc::ADC)
        return (adc.N == 1) ? [adc.T/2] .+ adc.delay : [range(0, adc.T; length=adc.N)...] .+ adc.delay
    end

    function koma_read(filename)
        path = @__DIR__
        seq = read_seq("$(path)/test_files/pulseq_read_comparison/$(filename).seq")
        N = length(seq)
        T0 = KomaMRIBase.get_block_start_times(seq)
        return Dict(
            "rf" => Dict(
                    "blocks"     => [i for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
                    "times"      => [get_theo_t_aux(seq[i].RF[1]) .+ T0[i] for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
                    "amplitudes" => [get_theo_A_aux(seq[i].RF[1]) for i in 1:N if KomaMRIBase.is_RF_on(seq[i])]
                ),
            "gx" => Dict(
                    "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
                    "times"      => [get_theo_t_aux(seq.GR[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
                    "amplitudes" => [get_theo_A_aux(seq.GR[1,i]) for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])]
                ),
            "gy" => Dict(
                    "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
                    "times"      => [get_theo_t_aux(seq.GR[2,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
                    "amplitudes" => [get_theo_A_aux(seq.GR[2,i]) for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])]
                ),
            "gz" => Dict(
                    "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
                    "times"      => [get_theo_t_aux(seq.GR[3,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
                    "amplitudes" => [get_theo_A_aux(seq.GR[3,i]) for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])]
                ),
            "adc" => Dict(
                    "blocks"     => [i for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
                    "times"      => [get_theo_t_aux(seq.ADC[i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
                )
        )
    end

    function pulseq_read(filename)
        path = @__DIR__
        seq_pulseq = matread("$(path)/test_files/pulseq_read_comparison/$(filename).mat")
        return Dict(
            "rf" => Dict(
                    "blocks"     => Int.(vec(seq_pulseq["rfBlocks"])),
                    "times"      => vec([vec(blockRfTimes) for blockRfTimes in seq_pulseq["rfTimes"]]),
                    "amplitudes" => vec([vec(blockRfAmplitudes) for blockRfAmplitudes in seq_pulseq["rfAmplitudes"]])
                ),
            "gx" => Dict(
                    "blocks"     => Int.(vec(seq_pulseq["gxBlocks"])),
                    "times"      => vec([vec(blockGxTimes) for blockGxTimes in seq_pulseq["gxTimes"]]),
                    "amplitudes" => vec([vec(blockGxAmplitudes) for blockGxAmplitudes in seq_pulseq["gxAmplitudes"]])
                ),
            "gy" => Dict(
                    "blocks"     => Int.(vec(seq_pulseq["gyBlocks"])),
                    "times"      => vec([vec(blockGyTimes) for blockGyTimes in seq_pulseq["gyTimes"]]),
                    "amplitudes" => vec([vec(blockGyAmplitudes) for blockGyAmplitudes in seq_pulseq["gyAmplitudes"]])
                ),
            "gz" => Dict(
                    "blocks"     => Int.(vec(seq_pulseq["gzBlocks"])),
                    "times"      => vec([vec(blockGzTimes) for blockGzTimes in seq_pulseq["gzTimes"]]),
                    "amplitudes" => vec([vec(blockGzAmplitudes) for blockGzAmplitudes in seq_pulseq["gzAmplitudes"]])
                ),
            "adc" => Dict(
                    "blocks"     => Int.(vec(seq_pulseq["adcBlocks"])),
                    "times"      => vec([vec(blockAdcTimes) for blockAdcTimes in seq_pulseq["adcTimes"]]),
                )
        )
    end

    function read_comparison(filename)
        # Define some names of the sequence dictionaries
        event_names = ["rf", "gx", "gy", "gz", "adc"]
        sample_names = ["times", "amplitudes"]
        # Read data for Koma and Pulseq
        seq_koma = koma_read(filename)      # it reads the .seq
        seq_pulseq = pulseq_read(filename)  # it reads the .mat (ground truth)
        # Checks for the blocks IDs of each event
        for event in event_names
            blocks_koma = seq_koma[event]["blocks"]
            blocks_pulseq = seq_pulseq[event]["blocks"]
            same_number_of_blocks = length(blocks_koma) == length(blocks_pulseq)
            @test same_number_of_blocks
            if same_number_of_blocks
                same_blocks = all(blocks_koma .== blocks_pulseq)
                @test same_blocks
                if !same_blocks
                    println("Mismatch for $(event) blocks")
                    data = [blocks_koma blocks_pulseq]
                    pretty_table(data; header=["koma", "pulseq"])
                end
            else
                println("Length mismatch for $(event):")
                println("Koma $(length(seq_koma[event])) != Pulseq $(length(seq_pulseq[event]))")
            end
        end
        # Checks for the timings of each event
        for event in event_names
            for sample in sample_names
                if !(event == "adc" && sample == "amplitudes")
                    for i in 1:length(seq_pulseq[event]["blocks"]) #1:5
                        samples_koma = seq_koma[event][sample][i]
                        samples_pulseq = seq_pulseq[event][sample][i]
                        same_points = all(.≈(abs.(samples_koma .- samples_pulseq), 0, atol=TOLERANCE))
                        @test same_points
                        if !same_points
                            println("Mismatch for $(event) $(sample). Block $(seq_pulseq[event]["blocks"][i])")
                            data = [samples_koma samples_pulseq]
                            pretty_table(data; header=["koma", "pulseq"])
                            # This is just an example about how to debug locally when there are problems
                            #if event == "gx" && sample == "amplitudes" #&& seq_pulseq[event]["blocks"][i] == 3
                            #    pretty_table(data; header=["koma", "pulseq"])
                            #    pretty_table([seq_koma[event]["times"][i] seq_pulseq[event]["times"][i]]; header=["koma_t", "pulseq_t"])
                            #    #pretty_table(data[end-10:end,:]; header=["koma", "pulseq"])
                            #end
                        end
                    end
                end
            end
        end
    end

    @testset "Trapezoidal_Grad" read_comparison("gr-trapezoidal")
    @testset "Uniformly_Shaped_Grad" read_comparison("gr-uniformly-shaped")
    @testset "Time_Shaped_Grad" read_comparison("gr-time-shaped")
    @testset "Pulse_RF" read_comparison("rf-pulse")
    @testset "Uniformly_Shaped_RF" read_comparison("rf-uniformly-shaped")
    @testset "FID" read_comparison("fid")
    @testset "Spiral" read_comparison("spiral")
    @testset "GRE" read_comparison("gre")
    @testset "EPI" read_comparison("epi")

    #@testset "cine_gre" read_comparison("cine_gre")                # success
    #@testset "DEMO_gre0" read_comparison("DEMO_gre0")              # success
    #@testset "DEMO_grep0" begin read_comparison("DEMO_grep0")      # success
    #@testset "epi" read_comparison("epi")                          # success
    #@testset "epi_lbl" read_comparison("epi_lbl")                  # success (takes long time)
    #@testset "epi_rs" read_comparison("epi_rs")                    # fail (cannot broadcast to a common size)
    #@testset "epi_rs_label" read_comparison("epi_rs_label")        # fail (cannot broadcast to a common size)
    #@testset "epi_se" read_comparison("epi_se")                    # success
    #@testset "epidiff_rs" read_comparison("epidiff_rs")            # fail (gr first)
    #@testset "epise_rs" read_comparison("epise_rs")                # fail (gr first)
    #@testset "fid" read_comparison("fid")                          # success
    #@testset "gre" read_comparison("gre")                          # success
    #@testset "gre3d" read_comparison("gre3d")                      # ??? (takes long time)
    #@testset "gre_gt" read_comparison("gre_gt")                    # success
    #@testset "gre_lbl" read_comparison("gre_lbl")                  # success
    #@testset "gre_rad" read_comparison("gre_rad")                  # success
    #@testset "gr-time-shaped" read_comparison("gr-time-shaped")    # success
    #@testset "gr-trapezoidal" read_comparison("gr-trapezoidal")    # success
    #@testset "MSE_test_KomaMRI" read_comparison("MSE_test_KomaMRI") # fail (cannot broadcast to a common size)
    #@testset "press" read_comparison("press")                      # success
    #@testset "rf-pulse" read_comparison("rf-pulse")                # success
    #@testset "rf-time-shaped" read_comparison("rf-time-shaped")    # success
    #@testset "rf-uniformly-shaped" read_comparison("rf-uniformly-shaped")  # success
    #@testset "selectiveRf" read_comparison("selectiveRf")          # fail (cannot broadcast to a common size)
    #@testset "spiral" read_comparison("spiral")                    # success
    #@testset "trufi" read_comparison("trufi")                      # fail (gr first)
    #@testset "tse" read_comparison("tse")                          # fail (gr first)
    #@testset "ute" read_comparison("ute")                          # success
    #@testset "ute_rs" read_comparison("ute_rs")                    # fail (cannot broadcast to a common size)
    #@testset "zte_petra" read_comparison("zte_petra")              # ??? (takes long time)

end
