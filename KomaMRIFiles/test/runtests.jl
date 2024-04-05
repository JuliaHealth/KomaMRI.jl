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
@testitem "gr-uniformly-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gr-uniformly-shaped") end
@testitem "gr-time-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gr-time-shaped") end
@testitem "rf-pulse" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("rf-pulse") end
#@testitem "rf-uniformly-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("rf-uniformly-shaped") end
#@testitem "fid" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("fid") end
#@testitem "spiral" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("spiral") end
#@testitem "gre" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre") end
#@testitem "epi" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi") end

#@testitem "cine_gre" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("cine_gre") end
#@testitem "DEMO_gre0" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("DEMO_gre0") end
#@testitem "DEMO_grep0" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("DEMO_grep0") end
#@testitem "epi" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi") end
#@testitem "epi_lbl" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi_lbl") end
#@testitem "epi_rs" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi_rs") end
#@testitem "epi_rs_label" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi_rs_label") end
#@testitem "epi_se" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epi_se") end
#@testitem "epidiff_rs" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epidiff_rs") end
#@testitem "epise_rs" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("epise_rs") end
#@testitem "fid" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("fid") end
#@testitem "gre" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre") end
#@testitem "gre3d" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre3d") end
#@testitem "gre_gt" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre_gt") end
#@testitem "gre_lbl" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre_lbl") end
#@testitem "gre_rad" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gre_rad") end
#@testitem "gr-time-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gr-time-shaped") end
#@testitem "gr-trapezoidal" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("gr-trapezoidal") end
#@testitem "MSE_test_KomaMRI" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("MSE_test_KomaMRI") end
#@testitem "press" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("press") end
#@testitem "rf-pulse" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("rf-pulse") end
#@testitem "rf-time-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("rf-time-shaped") end
#@testitem "rf-uniformly-shaped" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("rf-uniformly-shaped") end
#@testitem "selectiveRf" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("selectiveRf") end
#@testitem "spiral" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("spiral") end
#@testitem "trufi" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("trufi") end
#@testitem "ute" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("ute") end
#@testitem "ute_rs" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("ute_rs") end
#@testitem "zte_petra" tags=[:files] begin include(joinpath(@__DIR__, "test_files", "utils.jl")), read_comparison("zte_petra") end
