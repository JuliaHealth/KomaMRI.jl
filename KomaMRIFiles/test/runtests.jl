using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:files in ti.tags) #verbose=true

@testitem "Files" tags=[:files] begin
    using Suppressor

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

    # Test ReadPulseq
    @testset "ReadPulseq" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        seq = @suppress read_seq(path*"/test_files/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        # Test Pulseq compression-decompression
        shape = ones(100)
        num_samples, compressed_data = KomaMRIFiles.compress_shape(shape)
        shape2 = KomaMRIFiles.decompress_shape(num_samples, compressed_data)
        @test shape == shape2
    end

    # Test WritePulseq
    @testset "WritePulseq" begin

        # These two are solved by changing the precission in the RF comparison (now: atol=1e-5)
        # Note that just the angle of the RF is the problem
        # epi_rs. rfs positions: 1, 2, 68, 69, 135, 136, 202, 203, 269, 270. wrong: 69
        # epise_rs. rfs positions: 1, 2, 4, 61, 62, 64, 121, 122, 124. wrong: 4, 64, 124

        path = @__DIR__
        test_folder = joinpath(@__DIR__, "test_files", "pulseq")

        # Test for some .seq files
        filenames = ["DEMO_gre", "DEMO_grep", "epi_se", "epi", "external", "gre_rad",
            "spiral", "tabletop_tse_pulseq", "cine_gre", "epi_label", "epi_rs", "epise_rs"]
        for seq_filename_head in filenames
            seq_filename_head = seq_filename_head
            seq_original_filename = seq_filename_head * ".seq"
            seq_written_filename = seq_filename_head * "_written.seq"
            seq_original_file = joinpath(test_folder, seq_original_filename)
            seq_written_file = joinpath(test_folder, seq_written_filename)
            seq_original = @suppress read_seq(seq_original_file)
            write_seq(seq_original, seq_written_file)
            seq_written = @suppress read_seq(seq_written_file)
            rm(seq_written_file; force=true)
            @test seq_original ≈ seq_written
        end

    end

    @testset "read_seq_via_blocks_as_int_array" begin

        # Get the .seq file path
        path = @__DIR__
        test_folder = joinpath(@__DIR__, "test_files", "pulseq")

        filenames = ["bSSFP_FA30deg_TE10ms_TR20ms_2D_(69x64)_pulseq", "DEMO_gre",
            "DEMO_grep", "epi_se", "epi", "external", "gre_rad", "spiral", "tabletop_tse_pulseq",
            "cine_gre", "epi_label", "epi_rs", "epise_rs"]

        # Auxiliar functions to display speed
        time_read_seq(f) = @time @suppress read_seq(f)
        time_read_seq_via_blocks_as_int_array(f) = @time @suppress read_seq_via_blocks_as_int_array(f)

        # Force precompilation
        time_read_seq(joinpath(test_folder, "spiral.seq"))
        time_read_seq_via_blocks_as_int_array(joinpath(test_folder, "spiral.seq"))

        # Compare for all pulseq test sequences
        for seq_filename_head in filenames
            seq_filename = seq_filename_head * ".seq"
            seq_file = joinpath(test_folder, seq_filename)
            println("File: $(seq_filename)")
            seq = time_read_seq(seq_file)
            seq_blocks_no_id_sort = time_read_seq_via_blocks_as_int_array(seq_file)
            println()
            @test seq == seq_blocks_no_id_sort
        end

    end
end
