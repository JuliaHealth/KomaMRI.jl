using TestItems, TestItemRunner, KomaMRIBase

@run_package_tests filter=t_start->!(:skipci in t_start.tags)&&(:files in t_start.tags) #verbose=true

@testitem "Files" tags=[:files] begin
    using Suppressor, KomaMRIBase

    # Test Pulseq
    @testset "Pulseq" begin
        pth = @__DIR__
        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.5/gre_rad.seq") #Pulseq v1.5.1
        @test seq.DEF["FileName"] == "gre_rad.seq"
        @test seq.DEF["PulseqVersion"] == v"1.5.1"
        @test seq.DEF["signature"][:hash] == "80eae81bb6b808f2cb4ed5d23885009b"
        
        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.4/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "67ebeffe6afdf0c393834101c14f3990"

        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.4/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "efc5eb7dbaa82aba627a31ff689c8649"

        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.2/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] == v"1.2.1"
        @test seq.DEF["signature"][:hash] == "f291a24409c3e8de01ddb93e124d9ff2"

        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.2/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] == v"1.2.1"
        @test seq.DEF["signature"][:hash] == "e827cfff4436b65a6341a4fa0f6deb07"

        # Test Pulseq compression-decompression
        shape = ones(100)
        num_samples, compressed_data = KomaMRIFiles.compress_shape(shape)
        shape2 = KomaMRIFiles.decompress_shape(num_samples, compressed_data)
        @test shape == shape2

        # Test label capability
        using KomaMRIBase
        
        seq = @suppress read_seq(pth*"/test_files/pulseq/v1.4/label_test.seq") 
        label = get_label(seq)
        m = maximum(label)
        a = AdcLabels(4,0,0,0,0,0,0,2,0,0,0,0)
        bool = true
        for field in fieldnames(typeof(m))
            if getfield(m,field) != getfield(a,field)
                bool = false
                println(field)
            end
        end
        @test bool
    end
    # Test JEMRIS
    @testset "JEMRIS" begin
        pth = @__DIR__
        obj = read_phantom_jemris(pth*"/test_files/phantom/column1d.h5")
        @test obj.name == "column1d.h5"
    end
    # Test MRiLab
    @testset "MRiLab" begin
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_mrilab.mat"
        FRange_filename = pth * "/test_files/phantom/FRange.mat" #Slab within slice thickness
        obj = read_phantom_MRiLab(filename; FRange_filename)
        @test obj.name == "brain_mrilab.mat"
    end
    # Test Phantom (.phantom)
    @testset "Phantom" begin
        using KomaMRIBase
        pth = @__DIR__
        # NoMotion
        filename = pth * "/test_files/phantom/brain_nomotion_w.phantom"
        obj1 = brain_phantom2D()
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
    end
    @testset "SimpleAction" begin
        # SimpleAction
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_simplemotion_w.phantom"
        obj1 = brain_phantom2D()
        obj1.motion = MotionList(
            rotate(0.0, 0.0, 45.0, Periodic(period=1.0)),
            translate(0.0, 0.02, 0.0, TimeRange(t_start=0.0, t_end=0.5))
        )
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
    end
    @testset "ArbitraryAction" begin
        # ArbitraryAction
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_arbitrarymotion_w.phantom"
        obj1 = brain_phantom2D()
        Ns = length(obj1)
        K = 10
        t_start = 0.0
        t_end = 1.0
        obj1.motion = path(0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1), TimeRange(t_start, t_end)) 
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
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
    pth          = joinpath(@__DIR__, "test_files/pulseq/pulseq_read_comparison")
    pulseq_files = filter(endswith(".seq"), readdir(pth)) .|> x -> splitext(x)[1]
    for pulseq_file in pulseq_files
        #@show pulseq_file
        seq_koma   = @suppress read_seq("$pth/$pulseq_file.seq")
        seq_pulseq = matread("$pth/$pulseq_file.mat")["sequence"] .|> namedtuple
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
