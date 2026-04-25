using TestItems, TestItemRunner, KomaMRIBase

@run_package_tests filter=t_start->!(:skipci in t_start.tags)&&(:files in t_start.tags) #verbose=true

@testitem "JEMRIS" tags=[:files] begin
    using KomaMRIBase
    pth = @__DIR__
    obj = read_phantom_jemris(pth*"/test_files/phantom/column1d.h5")
    @test obj.name == "column1d.h5"
end

@testitem "MRiLab" tags=[:files] begin
    using KomaMRIBase
    pth = @__DIR__
    filename = pth * "/test_files/phantom/brain_mrilab.mat"
    FRange_filename = pth * "/test_files/phantom/FRange.mat" #Slab within slice thickness
    obj = read_phantom_MRiLab(filename; FRange_filename)
    @test obj.name == "brain_mrilab.mat"
end

@testitem "Phantom" tags=[:files] begin
    using KomaMRIBase
    @testset "NoMotion" begin
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_nomotion_w.phantom"
        obj1 = brain_phantom2D()
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
    end
    @testset "SimpleAction" begin
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_simplemotion_w.phantom"
        obj1 = brain_phantom2D()
        obj1.motion = MotionList(
            rotate(0.0, 0.0, 45.0, Periodic(period=1.0)),
            rotate(0.0, 0.0, 45.0, TimeRange(t_start=0.0, t_end=0.5), SpinRange(1:100); center=(0.0, 0.0, 1.0)),
            translate(0.0, 0.02, 0.0, TimeRange(t_start=0.0, t_end=0.5))
        )
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
    end
    @testset "ArbitraryAction" begin
        pth = @__DIR__
        filename = pth * "/test_files/phantom/brain_arbitrarymotion_w.phantom"
        obj1 = brain_phantom2D()
        Ns = length(obj1)
        K = 10
        t_start = 0.0
        t_end = 1.0
        obj1.motion = MotionList(
            path(    0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1),                      TimeRange(t_start, t_end)),
            flowpath(0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1), 0.01.*rand(Ns, K-1), rand(Bool, Ns, K-1), TimeRange(t_start, t_end))
        )
        write_phantom(obj1, filename)
        obj2 = read_phantom(filename)
        @test obj1 == obj2
    end
end

@testitem "Pulseq" tags=[:files, :pulseq] begin
    using MAT, KomaMRIBase, Suppressor
    function exported_sequence(data)
        buffer = IOBuffer()
        KomaMRIFiles.emit_pulseq(buffer, data)
        return KomaMRIBase.Sequence(KomaMRIFiles.read_seq_data(IOBuffer(take!(buffer))))
    end
    exported_sequence(seq::Sequence) = exported_sequence(@suppress KomaMRIFiles.write_seq_data(seq; check_timing=false))

    @testset "Basic Tests" begin
        pth = (@__DIR__ )*"/test_files/pulseq/basic_tests/"
        seq = @suppress read_seq(pth*"v1.5/gre_rad.seq") #Pulseq v1.5.1
        @test seq.DEF["FileName"] == "gre_rad.seq"
        @test seq.DEF["PulseqVersion"] == v"1.5.1"
        @test seq.DEF["signature"][:hash] == "80eae81bb6b808f2cb4ed5d23885009b"

        seq = @suppress read_seq(pth*"v1.5/unknown_ext.seq") #Pulseq v1.5.0 with unknown extensions
        @test seq.DEF["FileName"] == "unknown_ext.seq" 
        @test seq.DEF["PulseqVersion"] == v"1.5.0"
        
        seq = @suppress read_seq(pth*"v1.4/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "67ebeffe6afdf0c393834101c14f3990"

        seq = @suppress read_seq(pth*"v1.4/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "efc5eb7dbaa82aba627a31ff689c8649"

        seq = @suppress read_seq(pth*"v1.2/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] == v"1.2.1"
        @test seq.DEF["signature"][:hash] == "f291a24409c3e8de01ddb93e124d9ff2"

        seq = @suppress read_seq(pth*"v1.2/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] == v"1.2.1"
        @test seq.DEF["signature"][:hash] == "e827cfff4436b65a6341a4fa0f6deb07"

        sys = Scanner(GR_Δt=20e-6)
        seq = Sequence(sys)
        seq.DEF["GradientRasterTime"] = 10e-6
        @test_logs (:warn, r"GradientRasterTime") begin
            raster = KomaMRIFiles.PulseqRaster(seq, sys)
            @test raster.GradientRasterTime == sys.GR_Δt
        end

        sys = Scanner(B0=3.0, B1=17e-6, Gmax=40e-3, Smax=170.0, ADC_Δt=1e-6, DUR_Δt=10e-6, GR_Δt=10e-6, RF_Δt=1e-6, RF_ring_down_T=100e-6, RF_dead_time_T=10e-6, ADC_dead_time_T=10e-6)
        seq = Sequence(sys)
        seq.DEF["Name"] = "hardware-metadata"
        @addblock seq += RF(1e-6, 10e-6)
        data = @suppress KomaMRIFiles.write_seq_data(seq)
        buffer = IOBuffer()
        KomaMRIFiles.emit_pulseq(buffer, data)
        pulseq_text = String(take!(buffer))
        @test occursin("\nBlockDurationRaster ", pulseq_text)
        @test occursin("\nGradientRasterTime ", pulseq_text)
        for key in KomaMRIBase.PULSEQ_HW_DEFINITION_KEYS
            @test !occursin("\n$key ", pulseq_text)
        end

    end
    @testset "Compression-Decompression" begin
        shape = ones(100)
        num_samples, compressed_data = KomaMRIFiles.compress_shape(shape)
        shape2 = KomaMRIFiles.decompress_shape(num_samples, compressed_data)
        @test shape == shape2
    end
    @testset "Labels" begin
        pth = @__DIR__
        seq = @suppress read_seq(pth*"/test_files/pulseq/basic_tests/v1.4/label_test.seq") 
        label = get_labels(seq)
        m = maximum(label)
        a = AdcLabels(4,0,0,0,0,0,0,2,0,0,0,0)
        bool = true
        for field in fieldnames(typeof(m))
            if getfield(m,field) != getfield(a,field)
                bool = false
            end
        end
        @test bool
    end
    @testset "Read Comparison" begin
        namedtuple(x) = x[:]
        namedtuple(d::Dict) = (; (Symbol(k == "df" ? "Δf" : k) => namedtuple(v) for (k,v) in d)...)
        not_empty = ((ek, ep),) -> !isempty(ep.t)
        # Reading files
        pth          = joinpath(@__DIR__, "test_files/pulseq/read_comparison/")
        versions     = ["v1.2", "v1.3", "v1.4", "v1.5"]
        for v in versions
            pulseq_files = filter(endswith(".seq"), readdir(pth*v)) .|> x -> splitext(x)[1]
            for pulseq_file in pulseq_files
                #@show pulseq_file
                seq_koma   = @suppress read_seq("$pth$v/$pulseq_file.seq")
                seq_pulseq = matread("$pth$v/$pulseq_file.mat")["sequence"] .|> namedtuple
                @testset "$v/$pulseq_file" begin
                    for i in 1:length(seq_koma)
                        blk_koma   = get_samples(seq_koma, i)
                        blk_pulseq = NamedTuple{keys(blk_koma)}(seq_pulseq[i]) # Reorder keys
                        for (ev_koma, ev_pulseq) in Iterators.filter(not_empty, zip(blk_koma, blk_pulseq))
                            @test ev_koma.t ≈ ev_pulseq.t
                            @test ev_koma.A ≈ ev_pulseq.A
                        end
                    end
                end
            end
        end
    end
    @testset "Read Comparison Roundtrip" begin
        generated_prefix = "koma-generated-"
        pth = joinpath(@__DIR__, "test_files/pulseq/read_comparison/")
        versions = ["v1.2", "v1.3", "v1.4", "v1.5"]
        constraint_sys = Scanner(B1=1.0, Gmax=100.0, Smax=1e12, ADC_Δt=100e-9)
        mktempdir() do tmpdir
            for v in versions
                pulseq_files = filter(endswith(".seq"), readdir(pth * v)) .|> x -> splitext(x)[1]
                for pulseq_file in pulseq_files
                    label = "$v/$pulseq_file"
                    seq = @suppress read_seq("$pth$v/$pulseq_file.seq")
                    filename = joinpath(tmpdir, "$(generated_prefix)$(v)-$(pulseq_file).seq")
                    if occursin("JEMRIS", label) || occursin("gammaSTAR", label)
                        @suppress write_seq(seq, filename; check_timing=false)
                    else
                        @suppress write_seq(seq, filename; sys=constraint_sys)
                    end
                    seq2 = @suppress read_seq(filename)
                    @testset "$label" begin
                        @test seq2 ≈ seq
                    end
                end
            end
        end
    end
    @testset "Pulseq Event Representation Roundtrip" begin
        mktempdir() do tmpdir
            filename = joinpath(tmpdir, "trap.seq")
            seq = Sequence()
            seq += Grad(1e-3, 1e-3, 1e-5)
            @suppress write_seq(seq, filename)
            seq2 = @suppress read_seq(filename)
            @test seq2.GR[1, 1] isa KomaMRIBase.TrapezoidalGrad

            filename = joinpath(tmpdir, "uniform-rf.seq")
            Δt_rf = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
            seq = Sequence()
            seq += RF(ComplexF64[0.0, 1e-6, 0.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
            @suppress write_seq(seq, filename)
            seq2 = @suppress read_seq(filename)
            @test seq2.RF[1, 1] isa KomaMRIBase.UniformlySampledRF
        end
    end
    @testset "PulseqSequenceData Roundtrip" begin
        mktempdir() do tmpdir
            files = ("v1.4/rf-time-shaped.seq", "v1.5/gre.seq")
            for file in files
                src = joinpath(@__DIR__, "test_files/pulseq/read_comparison", file)
                data = @suppress KomaMRIFiles.read_seq_data(src)
                filename = joinpath(tmpdir, replace(file, "/" => "-"))
                @suppress write_seq(data, filename)
                seq2 = @suppress read_seq(filename)
                @test seq2 ≈ KomaMRIBase.Sequence(data)
            end
        end
    end
    @testset "Write Comparison" begin
        include("utils.jl")
        generated_prefix = "koma-generated-"
        signature_algorithms = ["md5", "sha1", "sha256"]
        mktempdir() do pth
            for (i, seq) in enumerate(round_trip_sequences())
                algorithm = signature_algorithms[mod1(i, length(signature_algorithms))]
                filename = joinpath(pth, "$(generated_prefix)$(seq.DEF["Name"]).seq")
                data = @suppress KomaMRIFiles.write_seq_data(seq; check_timing=false)
                qseq = exported_sequence(data)
                @suppress write_seq(data, filename; signatureAlgorithm=algorithm)
                seq2 = @suppress read_seq(filename)
                @test seq2 ≈ qseq
            end
        end
    end
    @testset "RF Compact Timing" begin
        raster = KomaMRIFiles.DEFAULT_RASTER
        Δt_rf = raster.RadiofrequencyRasterTime
        rf_event(seq) = begin
            prepared = @suppress KomaMRIFiles.prepare_pulseq_write(seq, raster)
            blocks, event_libraries = KomaMRIFiles.collect_pulseq_assets(prepared, raster)
            event_libraries.rf_library[blocks[1].rf_id]
        end

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(seq)
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id == 0
        @test event.delay ≈ Δt_rf
        @test event.center ≈ center + Δt_rf / 2
        @test qseq.RF[1, 1].center ≈ center

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], fill(Δt_rf, 2), 0.0, 2Δt_rf)
        qseq = exported_sequence(seq)
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.delay ≈ 2Δt_rf
        @test event.center ≈ center

        center = 0.83Δt_rf
        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2, center)
        event = rf_event(seq)
        @test event.time_shape_id == 0
        @test event.center ≈ center + Δt_rf / 2

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 0.25, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(seq)
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id == -1
        @test event.delay ≈ Δt_rf
        @test event.center ≈ center + Δt_rf / 2
        @test qseq.RF[1, 1].center ≈ center

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], [0.75Δt_rf, 1.25Δt_rf], 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(seq)
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.center ≈ center
        @test qseq.RF[1, 1].center ≈ center

        center = 0.83Δt_rf
        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], [0.75Δt_rf, 1.25Δt_rf], 0.0, 3Δt_rf / 2, center)
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.center ≈ center

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, Δt_rf / 4)
        event = rf_event(seq)
        @test event.time_shape_id > 0
    end
    @testset "RF Center Phase Representation" begin
        rf = RF(exp.(1im .* ([0.1, 0.2, 0.7] .+ 0.3)), 2.0, 0.0, 0.5, 1.0)
        @test rf.ϕ ≈ 0.5
        @test rf.A[2] ≈ 1.0 + 0.0im
        @test angle(cis(rf.ϕ) * rf.A[2]) ≈ rf.ϕ
    end
    @testset "Unsupported RF Frequency Waveform Write" begin
        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime, [0.0, 1.0], 0.0)
        prepared = @suppress KomaMRIFiles.prepare_pulseq_write(seq, KomaMRIFiles.DEFAULT_RASTER)
        @test_throws ArgumentError KomaMRIFiles.collect_pulseq_assets(prepared, KomaMRIFiles.DEFAULT_RASTER)
    end
    @testset "Legacy RF Center Fallback" begin
        seq = @suppress read_seq(joinpath(@__DIR__, "test_files/pulseq/read_comparison/v1.2/epi_JEMRIS.seq"))
        @test all(rf -> !is_RF_on(rf) || !isnothing(rf.center), vec(seq.RF))
    end
    @testset "Gradient Trap Preservation" begin
        raster = KomaMRIFiles.DEFAULT_RASTER
        trap = Sequence([Grad(1e-3, 2e-3, 1e-3, 1e-3, 0.5e-3)])
        filename = tempname() * ".seq"
        @suppress write_seq(trap, filename)
        roundtrip = @suppress read_seq(filename)
        prepared = @suppress KomaMRIFiles.prepare_pulseq_write(roundtrip, raster)
        _, libs = @suppress KomaMRIFiles.collect_pulseq_assets(prepared, raster)
        @test count(g -> g isa KomaMRIFiles.PulseqTrapGradEvent && !iszero(g.amplitude), values(libs.grad_library)) == 1
        @test all(g -> g isa KomaMRIFiles.PulseqTrapGradEvent, values(libs.grad_library))
        @test isempty(libs.shape_library)
    end
    @testset "RF Phase Shape Dedup" begin
        raster = KomaMRIFiles.DEFAULT_RASTER
        A = ComplexF64[0.0, 1.0, 1im, -1.0, 0.0]
        seq = Sequence()
        seq += RF(1e-6 .* A, 4raster.RadiofrequencyRasterTime)
        seq += RF(1e-6 .* A .* exp(1im * π / 3), 4raster.RadiofrequencyRasterTime)
        prepared = @suppress KomaMRIFiles.prepare_pulseq_write(seq, raster)
        _, libs = @suppress KomaMRIFiles.collect_pulseq_assets(prepared, raster)
        phase_ids = sort!(unique(rf.phase_id for rf in values(libs.rf_library)))
        @test length(phase_ids) == 1
    end
    @testset "RF Magnitude Shape Dedup" begin
        seq = @suppress read_seq(joinpath(@__DIR__, "test_files/pulseq/read_comparison/v1.5/gre.seq"))
        raster = KomaMRIFiles.PulseqRaster(seq)
        prepared = @suppress KomaMRIFiles.prepare_pulseq_write(seq, raster)
        _, libs = @suppress KomaMRIFiles.collect_pulseq_assets(prepared, raster)
        @test length(unique(rf.mag_id for rf in values(libs.rf_library))) == 1
        @test length(libs.shape_library) == 2
    end
end
