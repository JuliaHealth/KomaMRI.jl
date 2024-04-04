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

@testitem "gr-trapezoidal" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gr-trapezoidal") end
#@testitem "gr-uniformly-shaped" tags=[:files] read_comparison("gr-uniformly-shaped")
#@testitem "gr-time-shaped" tags=[:files] read_comparison("gr-time-shaped")
#@testitem "rf-pulse" tags=[:files] read_comparison("rf-pulse")
#@testitem "rf-uniformly-shaped" tags=[:files] read_comparison("rf-uniformly-shaped")
#@testitem "fid" tags=[:files] read_comparison("fid")
#@testitem "spiral" tags=[:files] read_comparison("spiral")
#@testitem "gre" tags=[:files] read_comparison("gre")
#@testitem "epi" tags=[:files] read_comparison("epi")

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
