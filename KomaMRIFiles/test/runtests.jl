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

end

@testitem "Pulseq compat" tags=[:files, :pulseq] begin
    using MAT, KomaMRIBase, Suppressor

    # Aux functions
    inside(x) = x[2:end-1]
    namedtuple(x) = x[:]
    namedtuple(d::Dict) = (; (Symbol(k == "df" ? "Δf" : k) => namedtuple(v) for (k,v) in d)...)
    not_empty = ((ek, ep),) -> !isempty(ep.t)

    # Reading files
    path         = joinpath(@__DIR__, "test_files/pulseq_read_comparison")
    pulseq_files = filter(endswith(".seq"), readdir(path)) .|> x -> splitext(x)[1]
    for pulseq_file in pulseq_files
        #@show pulseq_file
        seq_koma   = @suppress read_seq("$path/$pulseq_file.seq")
        seq_pulseq = matread("$path/$pulseq_file.mat")["sequence"] .|> namedtuple
        @testset "$pulseq_file" begin
            for i in 1:length(seq_koma)
                blk_koma   = get_samples(seq_koma, i)
                blk_pulseq = NamedTuple{keys(blk_koma)}(seq_pulseq[i]) # Reorder keys
                for (ev_koma, ev_pulseq) in Iterators.filter(not_empty, zip(blk_koma, blk_pulseq))
                    @test ev_koma.t ≈ ev_pulseq.t
                    @test inside(ev_koma.A) ≈ inside(ev_pulseq.A)
                    @test first(ev_koma.A)  ≈ first(ev_pulseq.A) || ev_koma.t[2] ≈ ev_koma.t[1]
                    @test last(ev_koma.A)   ≈ last(ev_pulseq.A)
                end
            end
        end
    end
end

@testitem "WritePulseq" tags=[:files, :pulseq] begin

    # These two are solved by changing the precission in the RF comparison (now: atol=1e-5)
    # Note that just the angle of the RF is the problem
    # epi_rs. rfs positions: 1, 2, 68, 69, 135, 136, 202, 203, 269, 270. wrong: 69
    # epise_rs. rfs positions: 1, 2, 4, 61, 62, 64, 121, 122, 124. wrong: 4, 64, 124

    path = @__DIR__
    test_folder = joinpath(@__DIR__, "test_files", "pulseq")

    # Test for some .seq files
    filenames = [
        "DEMO_gre",
        "DEMO_grep",
        "epi_se",
        "epi",
        "external",
        "gre_rad",
        "spiral",
        "tabletop_tse_pulseq",
        "cine_gre",
        "epi_label",
        "epi_rs",
        "epise_rs",
    ]
    for seq_filename_head in filenames
        seq_filename_head = seq_filename_head
        seq_original_filename = seq_filename_head * ".seq"
        seq_written_filename = seq_filename_head * "_written.seq"
        seq_original_file = joinpath(test_folder, seq_original_filename)
        seq_written_file = joinpath(test_folder, seq_written_filename)
        seq_original = read_seq(seq_original_file)
        write_seq(seq_original, seq_written_file)
        seq_written = read_seq(seq_written_file)
        rm(seq_written_file; force=true)
        @test seq_original ≈ seq_written
    end

end
