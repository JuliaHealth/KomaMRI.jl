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
    using MAT, KomaMRIBase
    function exported_sequence(data)
        buffer = IOBuffer()
        KomaMRIFiles.emit_pulseq(buffer, data)
        return KomaMRIBase.Sequence(KomaMRIFiles.read_seq_data(IOBuffer(take!(buffer))))
    end
    function pulseq_data_after_write(seq, filename="test.seq"; kwargs...)
        return mktempdir() do tmpdir
            path = joinpath(tmpdir, filename)
            write_seq(seq, path; verbose=false, kwargs...)
            read_seq_data(path; verbose=false)
        end
    end
    function matlab_signature_payload(filename)
        bytes = read(filename)
        sig_start = first(findfirst(b"[SIGNATURE]", bytes))
        @test bytes[sig_start - 1] == UInt8('\n')
        return bytes[1:(sig_start - 2)]
    end

    @testset "Basic Tests" begin
        pth = (@__DIR__ )*"/test_files/pulseq/basic_tests/"
        seq = read_seq(pth*"v1.5/gre_rad.seq"; verbose=false) #Pulseq v1.5.1
        @test seq.DEF["FileName"] == "gre_rad.seq"
        @test seq.DEF["PulseqVersion"] == v"1.5.1"
        @test seq.DEF["signature"][:hash] == "80eae81bb6b808f2cb4ed5d23885009b"

        seq = read_seq(pth*"v1.5/unknown_ext.seq"; verbose=false) #Pulseq v1.5.0 with unknown extensions
        @test seq.DEF["FileName"] == "unknown_ext.seq" 
        @test seq.DEF["PulseqVersion"] == v"1.5.0"
        
        seq = read_seq(pth*"v1.4/epi.seq"; verbose=false) #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "67ebeffe6afdf0c393834101c14f3990"
        @test all(rf -> !KomaMRIBase.is_on(rf) || rf.use isa KomaMRIBase.Excitation, seq.RF)

        seq = read_seq(pth*"v1.4/spiral.seq"; verbose=false) #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] == v"1.4.0"
        @test seq.DEF["signature"][:hash] == "efc5eb7dbaa82aba627a31ff689c8649"

        seq = read_seq(pth*"v1.2/epi_JEMRIS.seq"; verbose=false) #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] == v"1.2.1"
        @test seq.DEF["signature"][:hash] == "f291a24409c3e8de01ddb93e124d9ff2"

        seq = read_seq(pth*"v1.2/radial_JEMRIS.seq"; verbose=false) #Pulseq v1.2.1
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

        rf_raster = 1e-6
        block_raster = 10e-6
        gradient_raster = block_raster
        adc_raster = rf_raster
        rf_dead_time = block_raster
        adc_dead_time = block_raster
        rf_ring_down_time = 10block_raster
        rf_duration = block_raster
        rf_amp = 1e-6
        sys = Scanner(B0=3.0, B1=17rf_amp, Gmax=40e-3, Smax=170.0, ADC_Δt=adc_raster, DUR_Δt=block_raster, GR_Δt=gradient_raster, RF_Δt=rf_raster, RF_ring_down_time=rf_ring_down_time, RF_dead_time=rf_dead_time, ADC_dead_time=adc_dead_time)
        bad_rf_deadtime = Sequence(sys)
        @addblock bad_rf_deadtime += (RF(rf_amp, rf_duration), Duration(rf_dead_time + rf_duration + rf_ring_down_time))
        @test_throws ErrorException write_seq_data(bad_rf_deadtime; check_hw_limits=false, verbose=false)
        bad_rf_ringdown = Sequence(sys)
        @addblock bad_rf_ringdown += (RF(rf_amp, rf_duration, 0.0, rf_dead_time), Duration(rf_dead_time + rf_duration + rf_ring_down_time - block_raster))
        @test_throws ErrorException write_seq_data(bad_rf_ringdown; check_hw_limits=false, verbose=false)
        bad_rf_amplitude = Sequence(sys)
        @addblock bad_rf_amplitude += (RF(2sys.B1, rf_duration, 0.0, rf_dead_time), Duration(rf_dead_time + rf_duration + rf_ring_down_time))
        @test_throws ErrorException write_seq_data(bad_rf_amplitude; sys, verbose=false)
        bad_adc_deadtime = Sequence(sys)
        adc_samples = 2
        adc_dwell = 2adc_raster
        @addblock bad_adc_deadtime += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_dead_time), Duration(4block_raster))
        @test_throws ErrorException write_seq_data(bad_adc_deadtime; check_hw_limits=false, verbose=false)
        bad_adc_dwell = Sequence(sys)
        bad_dwell = adc_raster / 2
        @addblock bad_adc_dwell += (ADC(adc_samples, (adc_samples - 1) * bad_dwell, bad_dwell / 2), Duration(block_raster))
        @test_throws ErrorException write_seq_data(bad_adc_dwell; check_hw_limits=false, verbose=false)

        # With timing checks disabled, the writer still quantizes ADC dwell to the Pulseq raster.
        data = write_seq_data(bad_adc_dwell; check_timing=false, check_hw_limits=false, verbose=false)
        @test only(values(data.libraries.adc_library)).dwell == adc_raster

        # A block with one gradient axis should not allocate empty library events for the others.
        single_axis = Sequence(Scanner(Gmax=Inf, Smax=Inf))
        @addblock single_axis += (x=Grad(1e-3, 1e-3))
        data = write_seq_data(single_axis; verbose=false)
        @test only(data.blocks).gx_id != 0
        @test only(data.blocks).gy_id == 0
        @test only(data.blocks).gz_id == 0
        @test length(data.libraries.grad_library) == 1

        @testset "Write rejects events longer than block duration" begin
            # Strict Pulseq write duration checks use event end times, including sample-edge offsets.
            block_raster = 50e-9
            rf_raster = 100e-9
            grad_raster = 100e-9
            adc_raster = 100e-9
            rf_half_raster = rf_raster / 2
            grad_half_raster = grad_raster / 2
            fit_sys = Scanner(B1=Inf, Gmax=Inf, Smax=Inf, DUR_Δt=block_raster, GR_Δt=grad_raster, RF_Δt=rf_raster, ADC_Δt=adc_raster, RF_dead_time=0.0, RF_ring_down_time=0.0, ADC_dead_time=0.0)
            function expect_write_duration_error(seq, duration) # Strict write must not stretch this duration.
                seq.DUR[1] = duration
                @test_throws ErrorException KomaMRIFiles.write_seq_data(seq; sys=fit_sys, check_hw_limits=false, verbose=false)
            end
            function expect_write_duration_ok(seq, duration)
                seq.DUR[1] = duration
                @test !isnothing(KomaMRIFiles.write_seq_data(seq; sys=fit_sys, check_hw_limits=false, verbose=false))
            end

            # Compact sampled RF and gradients extend half a sample past their Koma event duration.
            compact_rf = Sequence(fit_sys)
            @addblock compact_rf += RF([1e-6, 1e-6], rf_raster, 0.0, rf_half_raster)
            compact_rf_end = dur(compact_rf) + rf_half_raster
            @test dur(compact_rf) ≈ rf_raster + rf_half_raster
            expect_write_duration_error(compact_rf, dur(compact_rf))
            expect_write_duration_ok(compact_rf, compact_rf_end)

            # Time-shaped events must reserve every explicit interval, including the last one.
            time_shaped_rf = Sequence(fit_sys)
            rf_intervals = [rf_raster, rf_raster]
            @addblock time_shaped_rf += RF([1e-6, 0.5e-6], rf_intervals, 0.0, 0.0)
            expect_write_duration_error(time_shaped_rf, first(rf_intervals))
            expect_write_duration_ok(time_shaped_rf, sum(rf_intervals))

            compact_grad = Sequence(fit_sys)
            @addblock compact_grad += (x=Grad([0.2e-3, 0.4e-3], grad_raster, grad_half_raster, grad_half_raster, 0.0))
            compact_grad_end = grad_raster + 2 * grad_half_raster
            expect_write_duration_error(compact_grad, compact_grad_end - grad_half_raster)
            expect_write_duration_ok(compact_grad, compact_grad_end)

            time_shaped_grad = Sequence(fit_sys)
            grad_intervals = [grad_raster, grad_raster]
            @addblock time_shaped_grad += (x=Grad([0.2e-3, 0.4e-3], grad_intervals))
            expect_write_duration_error(time_shaped_grad, first(grad_intervals))
            expect_write_duration_ok(time_shaped_grad, sum(grad_intervals))

            # ADC write duration includes the trailing half dwell after the last Koma sample.
            adc_samples = 2
            adc_dwell = adc_raster
            adc = ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_dwell / 2)
            adc_end = dur(adc) + adc_dwell / 2
            @test dur(adc) ≈ adc.delay + adc.T
            adc_edge = Sequence(fit_sys)
            @addblock adc_edge += adc
            expect_write_duration_error(adc_edge, dur(adc))
            expect_write_duration_ok(adc_edge, adc_end)
        end

    end
    @testset "Pulseq rotation extension" begin
        @testset "Unsupported required extensions warn" begin
            @test_logs (:warn, r"Ignoring unsupported required extension") KomaMRIFiles.read_extensions(
                IOBuffer("1 2 3\n\n"),
                "UNKNOWN",
                nothing,
                1,
                Dict{Int,Type{<:Extension}}(),
                Dict{Int,Dict{Int,Extension}}(),
                ["UNKNOWN"],
            )
        end

        mktempdir() do tmpdir
            @testset "ROTATIONS read/write roundtrip" begin
                filename = joinpath(tmpdir, "rotation.seq")
                q = QuaternionRot(cos(π / 4), 0, 0, sin(π / 4))
                seq = Sequence()
                @addblock seq += (
                    q,
                    x=Grad(1e-3, 1e-3, 0.1e-3, 0.1e-3),
                    y=Grad([0.2e-3, 0.4e-3, -0.1e-3], 0.9e-3, 0.0, 0.0, 0.1e-3),
                )
                write_seq(seq, filename; check_timing=false, verbose=false)
                pulseq_text = read(filename, String)
                @test occursin("RequiredExtensions ROTATIONS", pulseq_text)
                @test occursin("extension ROTATIONS", pulseq_text)

                raw = read_seq(filename; apply_rotations=false, verbose=false)
                @test any(ext -> ext isa QuaternionRot, raw.EXT[1])
                applied = read_seq(filename; verbose=false)
                @test any(ext -> ext isa QuaternionRot, applied.EXT[1])
                @test applied.DEF["RequiredExtensions"] == ["ROTATIONS"]

                encoded = apply_rotations(seq; reverse=true)
                t = collect(range(0, dur(seq); length=129))
                for (actual, expected) in zip(KomaMRIBase.get_grads(raw, t), KomaMRIBase.get_grads(encoded, t))
                    @test actual ≈ expected
                end
                for (actual, expected) in zip(KomaMRIBase.get_grads(applied, t), KomaMRIBase.get_grads(seq, t))
                    @test actual ≈ expected
                end

                rewritten = joinpath(tmpdir, "rotation_rewrite.seq")
                write_seq(applied, rewritten; check_timing=false, verbose=false)
                rewritten_text = read(rewritten, String)
                @test occursin("RequiredExtensions ROTATIONS", rewritten_text)
                @test occursin("extension ROTATIONS", rewritten_text)
            end

            @testset "Constant time-shaped gradients preserve area" begin
                filename = joinpath(tmpdir, "constant_time_shaped.seq")
                g = Grad(
                    [-0.2e-3, -0.2e-3, nextfloat(-0.2e-3)],
                    [1.19e-3, 2e-14], 0.25e-3, 0.25e-3,
                )
                seq = Sequence()
                @addblock seq += (y=g)
                write_seq(seq, filename; check_timing=false, verbose=false)
                roundtrip = read_seq(filename; verbose=false)
                @test area(roundtrip.GR[2,1]) ≈ area(g)
            end

            @testset "ROTATIONS gradient-library dedup" begin
                repeated = Sequence()
                gx = Grad([0.0, 0.2e-3, 0.4e-3, 0.2e-3], 12e-6, 2e-6, 2e-6)
                gy = Grad([0.0, -0.1e-3, -0.3e-3, -0.1e-3], 12e-6, 2e-6, 2e-6)
                gz = Grad([0.0, 0.3e-3, 0.1e-3, 0.3e-3], 12e-6, 2e-6, 2e-6)
                for rot in (QuaternionRot(rotz(π / 7)), QuaternionRot(rotx(-π / 5)), QuaternionRot(roty(π / 3) * rotz(π / 9)))
                    @addblock repeated += (rot, x=gx, y=gy, z=gz)
                end
                repeated_data = write_seq_data(apply_rotations(repeated); check_timing=false, check_hw_limits=false, verbose=false)
                @test length(repeated_data.libraries.grad_library) == 3
                @test all(g -> !(g isa KomaMRIFiles.PulseqArbGradEvent) || iszero(g.first) || abs(g.first) > abs(g.amplitude) * KomaMRIFiles.PULSEQ_SHAPE_ZERO_TOL, values(repeated_data.libraries.grad_library))
                @test all(g -> !(g isa KomaMRIFiles.PulseqArbGradEvent) || iszero(g.last) || abs(g.last) > abs(g.amplitude) * KomaMRIFiles.PULSEQ_SHAPE_ZERO_TOL, values(repeated_data.libraries.grad_library))
            end
        end

        @testset "MATLAB ROTATIONS fixture" begin
            fixture = joinpath(@__DIR__, "test_files/pulseq/basic_tests/v1.5/rotation_radial_tiny.seq")
            fixture_data = read_seq_data(fixture; verbose=false)
            @test KomaMRIFiles.supported_signature_digest(fixture_data.signature.type, matlab_signature_payload(fixture)) == fixture_data.signature.hash
            raw = read_seq(fixture; apply_rotations=false, verbose=false)
            applied = read_seq(fixture; verbose=false)
            @test applied.DEF["RequiredExtensions"] == ["ROTATIONS"]
            @test all(ext -> ext isa QuaternionRot, only.(raw.EXT))
            @test apply_rotations(raw) ≈ applied

            data = write_seq_data(applied; verbose=false)
            @test length(data.libraries.grad_library) <= 2
            @test count(g -> !iszero(g.amplitude), values(data.libraries.grad_library)) == 1
            @test length(data.libraries.adc_library) == 1
            @test length(data.libraries.extension_instance_library) == 3
            @test length(only(values(data.libraries.extension_spec_library))) == 3
        end
    end
    @testset "Compression-Decompression" begin
        shape = ones(100)
        num_samples, compressed_data = KomaMRIFiles.compress_shape(shape)
        shape2 = KomaMRIFiles.decompress_shape(num_samples, compressed_data)
        @test shape == shape2
    end
    @testset "Labels" begin
        mktempdir() do tmpdir
            # Extension chains must roundtrip without dropping labels from the same block.
            seq = Sequence()
            @addblock seq += (LabelSet(0, "SLC"), LabelSet(0, "LIN"))
            @addblock seq += (LabelSet(0, "SLC"), LabelSet(0, "LIN"), LabelInc(1, "SLC"))
            filename = joinpath(tmpdir, "label-chain.seq")
            write_seq(seq, filename; check_timing=false, verbose=false)
            roundtrip = read_seq(filename; verbose=false)
            @test roundtrip.EXT[1] == [LabelSet(0, "SLC"), LabelSet(0, "LIN")]
            @test roundtrip.EXT[2] == [LabelSet(0, "SLC"), LabelSet(0, "LIN"), LabelInc(1, "SLC")]
        end

        pth = @__DIR__
        seq = read_seq(pth*"/test_files/pulseq/basic_tests/v1.4/label_test.seq"; verbose=false)
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
                seq_koma   = read_seq("$pth$v/$pulseq_file.seq"; verbose=false)
                seq_pulseq = matread("$pth$v/$pulseq_file.mat")["sequence"] .|> namedtuple
                @testset "$v/$pulseq_file" begin
                    for i in 1:length(seq_koma)
                        blk_koma   = get_samples(seq_koma, i)
                        sample_keys = filter(k -> hasproperty(seq_pulseq[i], k), keys(blk_koma))
                        blk_koma = NamedTuple{sample_keys}(blk_koma)
                        blk_pulseq = NamedTuple{sample_keys}(seq_pulseq[i]) # Reorder keys
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
        mktempdir() do tmpdir
            for v in versions
                pulseq_files = filter(endswith(".seq"), readdir(pth * v)) .|> x -> splitext(x)[1]
                for pulseq_file in pulseq_files
                    label = "$v/$pulseq_file"
                    seq = read_seq("$pth$v/$pulseq_file.seq"; verbose=false)
                    filename = joinpath(tmpdir, "$(generated_prefix)$(v)-$(pulseq_file).seq")
                    # Legacy JEMRIS files omit AdcRasterTime but use 1 ns dwell precision.
                    occursin("JEMRIS", label) && (seq.DEF["AdcRasterTime"] = 1e-9)
                    write_seq(seq, filename; verbose=false)
                    seq2 = read_seq(filename; verbose=false)
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
            write_seq(seq, filename; verbose=false)
            seq2 = read_seq(filename; verbose=false)
            @test seq2.GR[1, 1] isa KomaMRIBase.TrapezoidalGrad

            @testset "Split-gradient block-edge continuity" begin
                filename = joinpath(tmpdir, "split-gradient.seq")
                seq = Sequence()
                gx1 = Grad(1e-3, 100e-6, 10e-6, 0.0, 0.0, 0.0, 1e-3)
                gx2 = Grad(1e-3, 100e-6, 0.0, 10e-6, 0.0, 1e-3, 0.0)
                @addblock seq += (x=gx1)
                @addblock seq += (x=gx2)
                write_seq(seq, filename; verbose=false)
                seq2 = read_seq(filename; verbose=false)
                @test seq2.GR[1, 1].last ≈ gx1.last
                @test seq2.GR[1, 2].first ≈ gx2.first
                t = collect(0:10e-6:dur(seq))
                @test KomaMRIBase.get_grads(seq2, t)[1] ≈ KomaMRIBase.get_grads(seq, t)[1]
                expected_slew = gx1.A / gx1.rise
                @test isnothing(check_hw_limits(seq2, Scanner(B1=Inf, Gmax=Inf, Smax=expected_slew, ADC_Δt=0.0)))
                @test_throws ErrorException check_hw_limits(seq2, Scanner(B1=Inf, Gmax=Inf, Smax=expected_slew / 2, ADC_Δt=0.0))
            end

            filename = joinpath(tmpdir, "uniform-rf.seq")
            Δt_rf = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
            seq = Sequence()
            seq += RF(ComplexF64[0.0, 1e-6, 0.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
            write_seq(seq, filename; verbose=false)
            seq2 = read_seq(filename; verbose=false)
            @test seq2.RF[1, 1] isa KomaMRIBase.UniformlySampledRF
        end
    end
    @testset "PulseqSequenceData Roundtrip" begin
        mktempdir() do tmpdir
            files = ("v1.4/rf-time-shaped.seq", "v1.5/gre.seq")
            for file in files
                src = joinpath(@__DIR__, "test_files/pulseq/read_comparison", file)
                data = KomaMRIFiles.read_seq_data(src; verbose=false)
                filename = joinpath(tmpdir, replace(file, "/" => "-"))
                write_seq(data, filename; verbose=false)
                seq2 = read_seq(filename; verbose=false)
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
                data = KomaMRIFiles.write_seq_data(seq; check_timing=false, verbose=false)
                qseq = exported_sequence(data)
                write_seq(data, filename; signatureAlgorithm=algorithm, verbose=false)
                seq2 = read_seq(filename; verbose=false)
                @test seq2 ≈ qseq
            end
        end
    end
    @testset "RF Compact Timing" begin
        Δt_rf = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
        rf_event(seq; check_timing=true) = begin
            data = pulseq_data_after_write(seq, "rf-compact.seq"; check_timing)
            data.libraries.rf_library[data.blocks[1].rf_id]
        end

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(KomaMRIFiles.write_seq_data(seq; verbose=false))
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id == 0
        @test event.delay ≈ Δt_rf
        @test event.center ≈ center + Δt_rf / 2
        @test qseq.RF[1, 1].center ≈ center

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], fill(Δt_rf, 2), 0.0, 2Δt_rf)
        qseq = exported_sequence(KomaMRIFiles.write_seq_data(seq; verbose=false))
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.delay ≈ Δt_rf
        @test event.center ≈ center + Δt_rf

        center = 0.83Δt_rf
        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2, center)
        event = rf_event(seq)
        @test event.time_shape_id == 0
        @test event.center ≈ center + Δt_rf / 2

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 0.25, 0.5, 1.0], 2Δt_rf, 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(KomaMRIFiles.write_seq_data(seq; verbose=false))
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.delay ≈ Δt_rf
        @test event.center ≈ center + Δt_rf / 2
        @test qseq.RF[1, 1].center ≈ center

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], [0.75Δt_rf, 1.25Δt_rf], 0.0, 3Δt_rf / 2)
        qseq = exported_sequence(KomaMRIFiles.write_seq_data(seq; verbose=false))
        center = qseq.RF[1, 1].center
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.center ≈ center + Δt_rf / 2
        @test qseq.RF[1, 1].center ≈ center

        center = 0.83Δt_rf
        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], [0.75Δt_rf, 1.25Δt_rf], 0.0, 3Δt_rf / 2, center)
        event = rf_event(seq)
        @test event.time_shape_id > 0
        @test event.center ≈ center + Δt_rf / 2

        seq = Sequence()
        seq += RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, 0.0, Δt_rf / 4)
        event = rf_event(seq; check_timing=false)
        @test event.time_shape_id > 0
    end
    @testset "RF Center Phase Representation" begin
        rf = RF(exp.(1im .* ([0.1, 0.2, 0.7] .+ 0.3)), 2.0, 0.0, 0.5, 1.0)
        @test rf.ϕ ≈ 0.5
        @test rf.A[2] ≈ 1.0 + 0.0im
        @test angle(cis(rf.ϕ) * rf.A[2]) ≈ rf.ϕ
    end
    @testset "Unsupported RF frequency waveform errors through write_seq" begin
        seq = Sequence()
        Δt_rf = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
        addblock!(seq, RF(ComplexF64[1.0, 0.5, 1.0], 2Δt_rf, [0.0, 1.0], 3Δt_rf / 2), Duration(10Δt_rf))
        mktempdir() do tmpdir
            @test_throws ArgumentError write_seq(seq, joinpath(tmpdir, "rf-frequency.seq"); verbose=false)
        end
    end
    @testset "Legacy RF Center Fallback" begin
        seq = read_seq(joinpath(@__DIR__, "test_files/pulseq/read_comparison/v1.2/epi_JEMRIS.seq"); verbose=false)
        @test all(rf -> !is_RF_on(rf) || !isnothing(rf.center), vec(seq.RF))
    end
    @testset "Zero-amplitude Pulseq RF" begin
        mktempdir() do dir
            filename = joinpath(dir, "zero-rf.seq")
            Δt_rf = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
            seq = Sequence()
            addblock!(seq, RF(ComplexF64[1.0, 0.5, 1.0], [Δt_rf, Δt_rf], 0.0, Δt_rf), Duration(10Δt_rf))
            write_seq(seq, filename; verbose=false)

            # Some Pulseq writers emit zero-amplitude RFs. Koma does not, so this test makes one manually.
            lines = readlines(filename)
            i = findfirst(line -> startswith(line, "1 ") && contains(line, " 1 2 3 "), lines)
            fields = split(lines[i])
            fields[2] = "0.0"
            lines[i] = join(fields, " ")
            write(filename, join(lines, "\n"))

            seq = read_seq(filename; verbose=false)
            @test !is_RF_on(seq.RF[1, 1])
        end
    end
    @testset "Trapezoidal gradient writes without shape libraries" begin
        trap = Sequence([Grad(1e-3, 2e-3, 1e-3, 1e-3, 0.5e-3)])
        libs = pulseq_data_after_write(trap, "trap.seq").libraries
        @test count(g -> g isa KomaMRIFiles.PulseqTrapGradEvent && !iszero(g.amplitude), values(libs.grad_library)) == 1
        @test all(g -> g isa KomaMRIFiles.PulseqTrapGradEvent, values(libs.grad_library))
        @test isempty(libs.shape_library)
    end

    @testset "RF event and shape dedup" begin
        rf_raster = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
        rf_duration = 2rf_raster
        rf_delay = 3rf_raster / 2
        rf_data(rfs) = begin
            seq = Sequence()
            for rf in rfs
                addblock!(seq, rf, Duration(20rf_raster))
            end
            pulseq_data_after_write(seq, "rf-dedup.seq")
        end
        rf_events(data) = [data.libraries.rf_library[block.rf_id] for block in data.blocks]

        rf_cases = (
            ("block", A -> RF(A[2], rf_duration, 0.0, rf_raster)),
            ("uniform", A -> RF(A, rf_duration, 0.0, rf_delay)),
            ("time-shaped", A -> RF(A, [rf_raster, rf_raster], 0.0, rf_raster)),
        )
        phases = [0.2, 0.6, 1.1]
        base_A = 1e-6 .* [1.0, 0.5, 0.25] .* cis.(phases)
        for (name, make_rf) in rf_cases
            @testset "$name RF phase/amplitude event dedup" begin
                rf = make_rf(base_A)
                data = rf_data((rf, 2rf, (-1) * rf, cis(π / 3) * rf, cis(2π) * rf))
                events = rf_events(data)
                @test data.blocks[1].rf_id == data.blocks[end].rf_id
                @test length(data.libraries.rf_library) == 4
                @test length(unique(rf.mag_id for rf in events)) == 1
                @test length(unique(rf.phase_id for rf in events)) == 1
                @test length(unique(rf.time_shape_id for rf in events)) == 1
            end
        end

        rf_time_cases = (
            ("block", (RF(1e-6, rf_duration, 0.0, rf_raster), RF(1e-6, 2rf_duration, 0.0, rf_raster))),
            ("uniform", (RF(base_A, rf_duration, 0.0, rf_delay), RF(base_A, 2rf_duration, 0.0, rf_duration))),
            ("time-shaped", (RF(base_A, [rf_raster, rf_raster], 0.0, rf_raster), RF(base_A, [rf_raster, rf_duration], 0.0, rf_raster))),
        )
        for (name, rfs) in rf_time_cases
            @testset "$name RF time shape split" begin
                events = rf_events(rf_data(rfs))
                @test length(events) == 2
                @test length(unique(rf.mag_id for rf in events)) == 1
                @test length(unique(rf.phase_id for rf in events)) == 1
                @test length(unique(rf.time_shape_id for rf in events)) == 2
            end
        end

        @testset "RF phase shape dedup across global phase" begin
            A1 = 1e-6 .* [1.0, 0.5, 0.25] .* cis.(phases)
            A2 = 1e-6 .* [1.0, 0.25, 0.5] .* cis.(phases .+ π / 3)
            events = rf_events(rf_data((RF(A1, rf_duration, 0.0, rf_delay), RF(A2, rf_duration, 0.0, rf_delay))))
            @test length(unique(rf.phase_id for rf in events)) == 1
            @test length(unique(rf.mag_id for rf in events)) == 2
        end

        @testset "RF magnitude shape dedup across phase shapes" begin
            magnitudes = 1e-6 .* [1.0, 0.5, 0.25]
            A1 = magnitudes .* cis.([0.1, 0.3, 0.7])
            A2 = magnitudes .* cis.([0.2, 0.9, 1.4])
            events = rf_events(rf_data((RF(A1, rf_duration, 0.0, rf_delay), RF(A2, rf_duration, 0.0, rf_delay))))
            @test length(unique(rf.mag_id for rf in events)) == 1
            @test length(unique(rf.phase_id for rf in events)) == 2
        end
    end

    @testset "Gradient event and shape dedup" begin
        grad_raster = KomaMRIFiles.DEFAULT_RASTER.GradientRasterTime
        grad_duration = 2grad_raster
        function write_gradient_libraries(grads)
            off = [Grad(0.0, 0.0) for _ in eachindex(grads)]
            GR = [reshape(grads, 1, :); reshape(off, 1, :); reshape(off, 1, :)]
            rf_events = reshape([RF(0.0, 0.0) for _ in eachindex(grads)], 1, :)
            ADCs = [ADC(0, 0.0) for _ in eachindex(grads)]
            seq = Sequence(GR, rf_events, ADCs, dur.(grads))
            return pulseq_data_after_write(seq, "gradient-dedup.seq").libraries
        end

        function test_gradient_dedup(name, grad)
            single = write_gradient_libraries([grad])
            repeated = write_gradient_libraries([scale * grad for scale in (1, -1, 2, -2, 1, -1, 2, -2)])
            @test length(repeated.grad_library) == 4
            @test length(repeated.shape_library) == length(single.shape_library)
            events = collect(values(repeated.grad_library))
            @test any(g -> g.amplitude < 0, events)
            @test any(g -> g.amplitude > 0, events)
            if name != "trap"
                @test length(unique(g.amp_shape_id for g in events)) == 1
                @test length(unique(g.time_shape_id for g in events)) == 1
            else
                @test all(g -> g isa KomaMRIFiles.PulseqTrapGradEvent, events)
            end
        end

        test_gradient_dedup("trap", Grad(0.2e-3, grad_duration, grad_raster, grad_raster, grad_raster))
        test_gradient_dedup("uniform", Grad([0.2e-3, 0.5e-3, -0.25e-3], grad_duration, 0.0, 0.0, grad_raster))
        test_gradient_dedup("time-shaped", Grad([0.2e-3, 0.5e-3, -0.25e-3], [grad_raster, grad_duration], 0.0, 0.0, grad_raster))

        @testset "uniform gradient time shape split" begin
            libs = write_gradient_libraries([
                Grad([0.2e-3, 0.5e-3, -0.25e-3], grad_duration, 0.0, 0.0, grad_raster),
                Grad([0.2e-3, 0.5e-3, -0.25e-3], 2grad_duration, 0.0, 0.0, grad_raster),
            ])
            events = collect(values(libs.grad_library))
            @test length(events) == 2
            @test length(unique(g.amp_shape_id for g in events)) == 1
            @test length(unique(g.time_shape_id for g in events)) == 2
        end

        @testset "time-shaped gradient time shape split" begin
            libs = write_gradient_libraries([
                Grad([0.2e-3, 0.5e-3, -0.25e-3], [grad_raster, grad_duration], 0.0, 0.0, grad_raster),
                Grad([0.2e-3, 0.5e-3, -0.25e-3], [grad_duration, grad_duration], 0.0, 0.0, grad_raster),
            ])
            events = collect(values(libs.grad_library))
            @test length(events) == 2
            @test length(unique(g.amp_shape_id for g in events)) == 1
            @test length(unique(g.time_shape_id for g in events)) == 2
        end
    end
    @testset "Phase-cycled RF and ADC event dedup" begin
        block_raster = KomaMRIFiles.DEFAULT_RASTER.BlockDurationRaster
        rf_raster = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
        adc_raster = KomaMRIFiles.DEFAULT_RASTER.AdcRasterTime
        adc_samples = 128
        adc_dwell = 160adc_raster
        seq = Sequence()
        A = fill(1e-6 + 0im, 3)
        rf_pulseq_delay = 140rf_raster
        rf = RF(A, 2rf_raster, 0.0, rf_pulseq_delay + rf_raster / 2)
        adc_delay = 480adc_raster
        adc = ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_delay)
        adc_end = delay(adc) - dwell(adc) / 2 + adc.N * dwell(adc)
        block_duration = ceil(adc_end / block_raster) * block_raster
        for i in 1:8
            phase = cis(π * (i - 1))
            addblock!(seq, phase * rf, phase * adc, Duration(block_duration))
        end
        libs = pulseq_data_after_write(seq, "phase-cycled-dedup.seq").libraries
        @test length(libs.rf_library) == 2
        @test length(libs.adc_library) == 2
    end
    @testset "bSSFP-style RF ADC dedup" begin
        block_raster = KomaMRIFiles.DEFAULT_RASTER.BlockDurationRaster
        rf_raster = KomaMRIFiles.DEFAULT_RASTER.RadiofrequencyRasterTime
        adc_raster = KomaMRIFiles.DEFAULT_RASTER.AdcRasterTime
        adc_samples = 128
        adc_dwell = 160adc_raster
        seq = Sequence()
        rf = RF(fill(1e-6 + 0im, 3), 2rf_raster, 0.0, block_raster + rf_raster / 2)
        adc = ADC(adc_samples, (adc_samples - 1) * adc_dwell, 180adc_raster)
        ramp_factors = range(1 / 15, stop=1, length=13)
        for (i, factor) in enumerate(ramp_factors)
            phase = cis(π * (i - 1))
            addblock!(seq, factor * phase * rf, Duration(16block_raster))
        end
        for i in 1:128
            phase = cis(π * (length(ramp_factors) + i - 1))
            addblock!(seq, phase * rf, phase * adc, Duration(206block_raster))
        end
        libs = pulseq_data_after_write(seq, "bssfp-dedup.seq").libraries
        @test length(libs.rf_library) == 14
        @test length(libs.adc_library) == 2
    end
end
