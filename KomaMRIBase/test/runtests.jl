using TestItems, TestItemRunner

@run_package_tests filter=t_start->!(:skipci in t_start.tags)&&(:base in t_start.tags) #verbose=true

@testitem "Sequence" tags=[:base] begin
    struct AddBlockTestOp end
    function Base.:*(::AddBlockTestOp, seq::Sequence)
        out = 3.0 * seq
        out.DEF["AddBlockTestOp"] = true
        return out
    end

    @testset "Init" begin
        sys = Scanner()
        B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
        EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0.0,0.0,0.0])
        @test dur(EX) ≈ durRF #RF length matches what is supposed to be

        #ACQ construction
        N = 101
        FOV = 23e-2
        EPI = PulseDesigner.EPI(FOV, N, sys)
        TE = 30e-3
        d1 = TE-dur(EPI)/2-dur(EX)
        d1 = d1 > 0 ? d1 : 0.0

        #Sequence construction
        seq = d1 > 0 ? EX + Delay(d1) + EPI : EX + EPI #Only add delay if d1 is positive (enough space)
        seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
        @test dur(seq) ≈ dur(EX) + d1 + dur(EPI) #Sequence duration matches what is supposed to be

        for nchannels in 1:5
            seq1 = Sequence(reshape([Grad(Float64(i), 1e-3) for i in 1:nchannels], nchannels, :))
            @test size(seq1.GR, 1) == max(3, nchannels)
        end
        seq1 = Sequence(reshape([Grad(1.0, 1e-3)], 1, :))
        @test seq1.GR[2,1] !== seq1.GR[3,1]
        seq1.GR[2,1].A = 2.0
        @test seq1.GR[3,1].A == 0.0
    end

    @testset "Sequence operations" begin
        seq1() = begin
            seq = Sequence(
                reshape([Grad(0.01, 0.01, 0.001)], 1, 1),
                reshape([RF(1e-6, 1e-3)], 1, 1),
                [ADC(4, 1e-3)],
            )
            seq.EXT[1] = [LabelSet(1, "LIN")]
            seq.DEF["nested"] = Dict("value" => 1)
            return seq
        end
        seq2() = Sequence(
            reshape([Grad(0.02, 0.02, 0.002)], 1, 1),
            reshape([RF(2e-6, 2e-3)], 1, 1),
            [ADC(8, 2e-3)],
        )

        left = seq1()
        right = seq2()
        seq = left + right
        @test seq !== left
        @test seq.GR.A ≈ [0.01 0.02; 0.0 0.0; 0.0 0.0]
        @test seq.RF.A ≈ [1e-6 2e-6]
        @test seq.ADC.N == [4; 8]
        @test (left - right).GR.A ≈ [left.GR.A -right.GR.A]

        seq.GR[1,1].A = 0.03
        seq.RF[1,1].A = 3e-6
        seq.ADC[1].N = 16
        seq.EXT[1][1] = LabelSet(2, "LIN")
        seq.DEF["nested"]["value"] = 2

        @test left.GR[1,1].A == 0.01
        @test left.RF[1,1].A ≈ 1e-6
        @test left.ADC[1].N == 4
        @test left.EXT[1][1] == LabelSet(1, "LIN")
        @test left.DEF["nested"]["value"] == 1

        seq.GR[1,2].A = 0.04
        seq.RF[1,2].A = 4e-6
        seq.ADC[2].N = 32
        @test right.GR[1,1].A == 0.02
        @test right.RF[1,1].A ≈ 2e-6
        @test right.ADC[1].N == 8

        left = seq1()
        right = seq2()
        batched = Sequence([left, Sequence(), right])
        @test batched.GR.A ≈ [left.GR.A right.GR.A]
        @test batched.RF.A ≈ [left.RF.A right.RF.A]
        @test batched.ADC.N == [left.ADC.N; right.ADC.N]
        batched.GR[1,1].A = 0.04
        batched.RF[1,1].A = 4e-6
        batched.ADC[1].N = 32
        @test left.GR[1,1].A == 0.01
        @test left.RF[1,1].A ≈ 1e-6
        @test left.ADC[1].N == 4

        block = Sequence([Grad(0.01, 0.01, 0.001)])
        empty = Sequence()
        from_empty = Sequence([empty, block])
        @test from_empty !== empty
        from_empty.GR[1,1].A = 0.02
        @test length(empty) == 0
        @test block.GR[1,1].A == 0.01

        left = seq1()
        right = seq2()
        appended = Sequence()
        @test append!(appended, left) === appended
        @test push!(appended, right) === appended
        @test appended.GR.A ≈ [left.GR.A right.GR.A]
        @test appended.RF.A ≈ [left.RF.A right.RF.A]
        @test appended.ADC.N == [left.ADC.N; right.ADC.N]
        appended.GR[1,1].A = 0.05
        appended.RF[1,1].A = 5e-6
        appended.ADC[1].N = 64
        appended.GR[1,2].A = 0.06
        @test left.GR[1,1].A == 0.01
        @test left.RF[1,1].A ≈ 1e-6
        @test left.ADC[1].N == 4
        @test right.GR[1,1].A == 0.02

        g = Grad(0.04, 0.01, 0.001)
        rf = RF(4e-6, 1e-3)
        adc = ADC(16, 1e-3)
        seqg_left = seq1()
        seqrf_right = seq1()
        seqadc_left = seq1()
        seqg = Sequence([seqg_left])
        seqrf = Sequence([seqrf_right])
        seqadc = Sequence([seqadc_left])
        addblock!(seqg; x=g)
        addblock!(seqrf, rf)
        addblock!(seqadc, adc)
        seqg.GR[1,end].A = 0.05
        seqrf.RF[1,end].A = 5e-6
        seqadc.ADC[end].N = 32
        @test g.A == 0.04
        @test seqrf_right.RF[1,1].A ≈ 1e-6
        @test adc.N == 16

        seq0 = Sequence(
            [Grad(0.01, 0.01, 0.001); Grad(0.02, 0.01, 0.001); Grad(0.03, 0.01, 0.001);;],
            reshape([RF(1e-6, 1e-3)], 1, 1),
            [ADC(4, 1e-3)],
        )
        tuple_block = seq0 + (rf, adc)
        @test tuple_block.RF[1,end].A ≈ rf.A
        @test tuple_block.ADC[end].N == adc.N
        optional_seq = Sequence()
        for acquired in (false, true)
            @addblock optional_seq += rotz(0.0) * (x=g, acquired ? adc : nothing)
        end
        @test optional_seq.ADC.N == [0, adc.N]
        @test_throws MethodError seq0 + (0, 1, 0)
        @test_throws ErrorException seq0 + (rf, g)

        scaled = 2.0 * seq0
        @test scaled !== seq0
        @test scaled.GR.A ≈ 2 .* seq0.GR.A
        @test scaled.RF[1,1].A ≈ seq0.RF[1,1].A
        @test scaled.ADC[1].ϕ == seq0.ADC[1].ϕ

        phased = im * seq0
        @test phased !== seq0
        @test phased.GR.A ≈ seq0.GR.A
        @test phased.RF[1,1].ϕ ≈ π / 2
        @test phased.ADC[1].ϕ ≈ π / 2

        rotated = rotz(π / 2) * seq0
        @test rotated !== seq0
        @test rotated.GR[1,1].A ≈ -seq0.GR[2,1].A
        @test rotated.GR[2,1].A ≈ seq0.GR[1,1].A
        @test rotated.GR[3,1].A ≈ seq0.GR[3,1].A
        @test rotated.RF[1,1].A ≈ seq0.RF[1,1].A
        @test rotated.ADC[1].N == seq0.ADC[1].N

        @testset "Sequence k-space rotation" begin
            off_grad = Grad(0.0, 0.0)
            delayed_trap = Grad(0.4e-3, 0.9e-3, 0.1e-3, 0.15e-3, 0.05e-3)
            edge_trap = Grad(0.3e-3, 0.7e-3, 0.0, 0.0, 0.0)
            delayed_uniform = Grad([0.1e-3, 0.5e-3, -0.2e-3], 0.9e-3, 0.0, 0.0, 0.17e-3)
            ramped_time = Grad([0.2e-3, -0.1e-3, 0.3e-3], [0.2e-3, 0.4e-3], 0.07e-3, 0.0, 0.0)
            edge_time = Grad([0.2e-3, -0.1e-3, 0.3e-3], [0.2e-3, 0.4e-3], 0.0, 0.0, 0.0, 0.1e-3, 0.0)
            rotation_cases = (
                (delayed_trap, off_grad, edge_trap),
                (edge_trap, delayed_uniform, off_grad),
                (delayed_uniform, ramped_time, delayed_trap),
                (ramped_time, edge_time, edge_trap),
                (edge_time, delayed_trap, delayed_uniform),
            )
            rotations = (
                [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
                rotx(π), roty(π), rotz(π),
                rotx(2π), roty(2π), rotz(2π),
                rotz(π / 4),
                rotx(π / 7) * roty(-π / 5) * rotz(π / 6),
            )
            for gradients in rotation_cases
                adc_duration = any(!iszero(gr.last) for gr in gradients) ?
                    maximum(dur, gradients) :
                    maximum(dur, gradients) + 0.2e-3
                seq = Sequence()
                addblock!(seq, ADC(257, adc_duration); x=gradients[1], y=gradients[2], z=gradients[3])
                _, kspace = get_kspace(seq)
                for R in rotations
                    _, rotated_kspace = get_kspace(R * seq)
                    @test rotated_kspace ≈ kspace * R'
                end
            end
        end

        negated = -seq0
        divided = seq0 / 2
        @test negated.GR.A ≈ -seq0.GR.A
        @test divided.GR.A ≈ seq0.GR.A ./ 2
        @test negated.RF[1,1].A ≈ seq0.RF[1,1].A
        @test divided.ADC[1].N == seq0.ADC[1].N
        @testset "Sequence transforms copy events" begin
            for transformed in (scaled, phased, rotated, negated, divided)
                @test transformed.GR[1,1] !== seq0.GR[1,1]
                @test transformed.RF[1,1] !== seq0.RF[1,1]
                @test transformed.ADC[1] !== seq0.ADC[1]
            end
        end

        @test_throws MethodError seq0 * 2.0
        @test_throws MethodError seq0 * im
        @test_throws MethodError seq0 * rotz(π / 2)

        for seq in (scaled, phased, negated, divided)
            seq.GR[1,1].A = 0.04
            seq.RF[1,1].A = 4e-6
            seq.ADC[1].N = 32
        end
        rotated.GR[1,1].delay += 1e-6
        rotated.RF[1,1].A = 4e-6
        rotated.ADC[1].N = 32
        @test seq0.GR[1,1].A == 0.01
        @test seq0.RF[1,1].A ≈ 1e-6
        @test seq0.ADC[1].N == 4

        g_trap = Grad(1.0, 1.0)
        g_uniform = Grad([1.0, 2.0], 1.0)
        g_time = Grad([1.0, 2.0], [1.0])
        rf_block = RF(1.0, 1.0)
        rf_uniform = RF(ComplexF64[1.0, 2.0], 1.0)
        rf_time = RF(ComplexF64[1.0, 2.0], [1.0])
        rf_freq = RF(ComplexF64[1.0, 2.0], 1.0, [0.0, 1.0])
        mixed = Sequence([
            Sequence(reshape([g], 1, 1), reshape([rf], 1, 1))
            for g in (g_trap, g_uniform, g_time)
            for rf in (rf_block, rf_uniform, rf_time, rf_freq)
        ])
        @test eltype(mixed.GR) == Union{typeof(g_trap),typeof(g_uniform),typeof(g_time)}
        @test eltype(mixed.RF) == Union{typeof(rf_block),typeof(rf_uniform),typeof(rf_time),typeof(rf_freq)}

        macroseq = Sequence()
        counter = 0
        @addblock begin
            macroseq += (rf, z=g)
            macroseq += seq2()
            counter += 1
        end
        @test counter == 1
        @test length(macroseq) == 2
        @test macroseq.RF[1,1].A ≈ rf.A
        @test macroseq.GR[3,1].A ≈ g.A
        @test macroseq.GR[1,2].A ≈ seq2().GR[1,1].A

        axis_kw = (x=g, z=2g)
        macroseq = Sequence()
        @addblock macroseq += (rf; axis_kw...) + (; axis_kw...)
        @test length(macroseq) == 2
        @test macroseq.RF[1,1].A ≈ rf.A
        @test macroseq.GR[1,1].A ≈ g.A
        @test macroseq.GR[3,1].A ≈ (2g).A
        @test macroseq.GR[1,2].A ≈ g.A
        @test macroseq.GR[3,2].A ≈ (2g).A

        gxseq = Sequence()
        @addblock gxseq += (x=g)
        @test length(gxseq) == 1
        @test gxseq.GR[1,1].A ≈ g.A

        mixedseq = Sequence()
        left_chunk = seq1()
        right_chunk = seq2()
        @addblocks begin
            mixedseq += left_chunk + (rf, x=g) + right_chunk
        end
        @test length(mixedseq) == 3
        @test mixedseq.GR[1,1].A ≈ left_chunk.GR[1,1].A
        @test mixedseq.RF[1,2].A ≈ rf.A
        @test mixedseq.GR[1,2].A ≈ g.A
        @test mixedseq.GR[1,3].A ≈ right_chunk.GR[1,1].A
        mixedseq.RF[1,2].A = 7e-6
        mixedseq.GR[1,2].A = 0.07
        @test rf.A ≈ 4e-6
        @test g.A == 0.04

        spliced = Sequence()
        contents = (LabelSet(0, "LIN"), LabelInc(1, "ECO"))
        @addblocks begin
            spliced += (contents...)
            spliced += (contents..., Delay(20e-3), x=g)
        end
        @test length(spliced) == 2
        @test spliced.EXT[1] == Extension[contents...]
        @test spliced.EXT[2] == Extension[contents...]
        @test spliced.GR[1,2].A ≈ g.A
        @test dur(spliced[2]) ≈ 20e-3

        bssp = Sequence()
        @addblock bssp += (rf, z=g)
        @test length(bssp) == 1
        @test bssp.RF[1,1].A ≈ rf.A
        @test bssp.GR[3,1].A ≈ g.A
        bssp.RF[1,1].A = 8e-6
        bssp.GR[3,1].A = 0.08
        @test rf.A ≈ 4e-6
        @test g.A == 0.04

        transformed_tuple = Sequence()
        @addblock transformed_tuple += rotz(π / 2) * (x=g, adc) + im * (rf, z=g) + 0.5 * (x=g, adc)
        @test length(transformed_tuple) == 3
        @test transformed_tuple.GR[1,1].A ≈ 0 atol=1e-15
        @test transformed_tuple.GR[2,1].A ≈ g.A
        @test transformed_tuple.ADC[1].N == adc.N
        @test transformed_tuple.RF[1,2].ϕ ≈ π / 2
        @test transformed_tuple.GR[3,2].A ≈ g.A
        @test transformed_tuple.GR[1,3].A ≈ 0.5g.A
        transformed_tuple.GR[2,1].A = 0.09
        transformed_tuple.RF[1,2].A = 9e-6
        transformed_tuple.ADC[3].N = 64
        @test g.A == 0.04
        @test rf.A ≈ 4e-6
        @test adc.N == 16

        transformed_block = Sequence()
        @addblock transformed_block += (x=g, adc)
        transformed_sequence = Sequence()
        @addblock transformed_sequence += AddBlockTestOp() * transformed_block
        transformed_tuple = Sequence()
        @addblock transformed_tuple += AddBlockTestOp() * (x=g, adc)
        @test transformed_tuple.GR[1,1].A ≈ transformed_sequence.GR[1,1].A
        @test transformed_tuple.ADC[1].N == transformed_sequence.ADC[1].N
        @test transformed_tuple.DEF["AddBlockTestOp"] == true
        transformed_tuple.GR[1,1].A = 0.12
        @test g.A == 0.04

        bssp = Sequence()
        line = 0
        @addblocks for _ in 1:3
            bssp += (rf, z=g)
            line += 1
        end
        @test length(bssp) == 3
        @test line == 3
        @test all(block -> block.RF[1,1].A ≈ rf.A, bssp)

        default_sys = Scanner()
        block_raster = default_sys.DUR_Δt
        rf_raster = default_sys.RF_Δt
        off_raster_duration = block_raster + rf_raster / 2

        checked = Sequence()
        @addblock check_timing=true check_hw_limits=true checked += (RF(1e-6, block_raster), x=Grad(1e-3, block_raster))
        @test length(checked) == 1
        checks = (; check_timing=true, check_hw_limits=true)
        @addblock checks checked += (RF(1e-6, block_raster), x=Grad(1e-3, block_raster))
        @test length(checked) == 2
        @test_throws ErrorException @addblock check_timing=true checked += RF(1e-6, off_raster_duration)
        @test length(checked) == 2
        @test_throws ErrorException @addblock check_timing=true sys=Scanner(DUR_Δt=2block_raster) checked += RF(1e-6, block_raster)
        @test length(checked) == 2

        hw_checked = Sequence()
        @addblock check_hw_limits=true hw_checked += RF(1e-6, off_raster_duration)
        @test length(hw_checked) == 1
        @test_throws ErrorException @addblock check_hw_limits=true sys=Scanner(B1=0.5e-6) hw_checked += RF(1e-6, block_raster)
        @test length(hw_checked) == 1

        checked_loop = Sequence()
        @addblocks check_timing=true check_hw_limits=true for _ in 1:2
            checked_loop += (RF(1e-6, block_raster), x=Grad(1e-3, block_raster))
        end
        @test length(checked_loop) == 2
        @addblocks checks for _ in 1:2
            checked_loop += (RF(1e-6, block_raster), x=Grad(1e-3, block_raster))
        end
        @test length(checked_loop) == 4
        @test_throws ErrorException @addblocks check_timing=true for _ in 1:1
            checked_loop += RF(1e-6, off_raster_duration)
        end
        @test length(checked_loop) == 4

        @test_throws ErrorException addblock!(Sequence(), rf, g)
        @test_throws ErrorException macroexpand(@__MODULE__, :(@addblock badseq = (RF(1e-6, $block_raster))))
    end

    @testset "Grad" begin
        #Testing gradient concatenation, breakes in some Julia versions
        A1, A2, T = rand(3)
        g1, g2 = Grad(A1,T), Grad(A2,T)
        GR = [g1;g2;;]
        GR2 = reshape([g1;g2],:,1)
        @test GR.A ≈ GR2.A

        g_uniform = Grad([A1, A2], T)
        g_time = Grad([A1, A2], [T])
        GR = [g1; g_uniform; g_time;;]
        @test size(GR) == (3, 1)
        @test GR[1] === g1
        @test GR[2] === g_uniform
        @test GR[3] === g_time
        @test eltype(Sequence(GR).GR) == Union{typeof(g1),typeof(g_uniform),typeof(g_time)}

        #Sanity checks of contructors (A [T], T[s], rise[s], fall[s], delay[s])
        A, T = 0.1, 1e-3
        grad = Grad(A, T)

        A, T = rand(2)
        g1, g2 = Grad(A,T), Grad(A,T,0.0,0.0,0.0)
        @test g1 ≈ g2

        A, T, ζ = rand(3)
        g1, g2 = Grad(A,T,ζ), Grad(A,T,ζ,ζ,0.0)
        @test g1 ≈ g2

        A, T, delay, ζ = rand(4)
        g1, g2 = Grad(A,T,ζ,delay), Grad(A,T,ζ,ζ,delay)
        @test g1 ≈ g2

        # Test construction with shape function
        T, N = 1e-3, 100
        f = t -> sin(π*t / T)
        gradw = Grad(f, T, N)
        @test gradw.A ≈ f.(range(0.0, T; length=N))

        # Test Grad operations
        α = 3
        gradt = α * grad
        @test gradt.A ≈ α * grad.A
        gradt = grad * α
        @test gradt.A ≈ α * grad.A
        gradt = grad / α
        @test gradt.A ≈ grad.A / α
        grads = grad + gradt
        @test grads.A ≈ grad.A + gradt.A
        A1, A2, A3 = 0.1, 0.2, 0.3
        v1 = [Grad(A1,T); Grad(A2,T); Grad(A3,T)]
        v2 = [Grad(A2,T); Grad(A3,T); Grad(A1,T)]
        v3 = v1 + v2
        @test [v3[i].A for i=1:length(v3)] ≈ [v1[i].A + v2[i].A for i=1:length(v1)]
        gradr = grad - gradt
        @test gradr.A ≈ grad.A - gradt.A

        @testset "Gradient addition preserves types for matching timing" begin
            same_timing_pairs = (
                (Grad(0.4, 0.8e-3, 0.1e-3, 0.2e-3, 0.05e-3), Grad(-0.1, 0.8e-3, 0.1e-3, 0.2e-3, 0.05e-3), KomaMRIBase.TrapezoidalGrad),
                (Grad([0.1, 0.5, -0.2], 0.9e-3, 0.1e-3, 0.15e-3, 0.02e-3), Grad([0.3, -0.2, 0.4], 0.9e-3, 0.1e-3, 0.15e-3, 0.02e-3), KomaMRIBase.UniformlySampledGrad),
                (Grad([0.1, 0.5, -0.2], [0.3e-3, 0.6e-3], 0.1e-3, 0.15e-3, 0.02e-3), Grad([0.3, -0.2, 0.4], [0.3e-3, 0.6e-3], 0.1e-3, 0.15e-3, 0.02e-3), KomaMRIBase.TimeShapedGrad),
            )
            for (ga, gb, grad_type) in same_timing_pairs
                gsum = ga + gb
                @test gsum isa grad_type
                @test gsum.A ≈ ga.A .+ gb.A
                @test gsum.T == ga.T
                @test gsum.rise == ga.rise
                @test gsum.fall == ga.fall
                @test gsum.delay == ga.delay
            end

            constant_sum = Grad(fill(1.0, 4), 1e-3) + Grad(2.0, 1e-3)
            @test constant_sum isa KomaMRIBase.TrapezoidalGrad
            @test constant_sum.A ≈ 3.0
            @test length(ampls(constant_sum)) == 4
        end

        @testset "Gradient addition matches sampled waveform for mixed timing" begin
            _sample_xgrad(g, t) = first(KomaMRIBase.get_grads(Sequence(reshape([g], 1, 1)), t))
            ga = Grad(0.4, 0.8e-3, 0.1e-3, 0.2e-3, 0.05e-3)
            gb = Grad([0.1, 0.5, -0.2], 0.9e-3, 0.12e-3, 0.08e-3, 0.17e-3)
            t = sort!(unique!(vcat(times(ga), times(gb))))
            @test isapprox(_sample_xgrad(ga + gb, t), _sample_xgrad(ga, t) .+ _sample_xgrad(gb, t); atol=1e-12)
            @test area(ga + gb) ≈ area(ga) + area(gb)
        end

        gradt = -grad
        @test gradt.A ≈ -grad.A
        vc = vcat(v1, v2)
        @test [vc[1,j].A for j=1:length(v1)] ≈ [v1[i].A for i=1:length(v1)]
        @test [vc[2,j].A for j=1:length(v2)] ≈ [v2[i].A for i=1:length(v2)]
        vc = vcat(v1, v2, v3)
        @test [vc[1,j].A for j=1:length(v1)] ≈ [v1[i].A for i=1:length(v1)]
        @test [vc[2,j].A for j=1:length(v2)] ≈ [v2[i].A for i=1:length(v2)]
        @test [vc[3,j].A for j=1:length(v3)] ≈ [v3[i].A for i=1:length(v3)]
        delay, rise, T, fall = 1e-6, 2e-6, 10e-3, 3e-6
        gr = Grad(A, T, rise, fall, delay)
        @test dur(gr) ≈ delay + rise + T + fall
        @test area(Grad(2.0, 3.0, 1.0, 1.0, 0.0)) ≈ 8.0
        @test area(Grad([0.0, 1.0, 0.0], 2.0)) ≈ 1.0
        @test area(Grad([0.0, 1.0, 0.0], 2.0, 0.5, 0.5, 0.0, -0.5, -0.5)) ≈ 0.75
        @test area(Grad([0.0, 1.0, 0.0], [1.0, 1.0])) ≈ 1.0
        T1, T2, T3 = 1e-3, 2e-3, 3e-3
        vt = [Grad(A1,T1); Grad(A2,T2); Grad(A3,T3)]
        @test dur(vt) ≈ [maximum([T1, T2, T3])]

        # Just checking to ensure that show() doesn't get stuck and that it is covered
        show(IOBuffer(), "text/plain", grad)
        @test true

    end

    @testset "RF" begin
        #Testing gradient concatenation, breakes in some Julia versions
        A1, A2, T = rand(3)
        r1, r2 = RF(A1,T), RF(A2,T)
        R = [r1;r2;;]
        R2 = reshape([r1;r2],:,1)
        @test R.A ≈ R2.A

        #Sanity checks of constructors (A [T], T [s], Δf[Hz], delay [s])
        A, T = rand(2)
        r1 = RF(A,T)
        r2 = RF(A,T,0.0,0.0,r1.center)
        @test r1 ≈ r2
        A, T, Δf = rand(3)
        r1 = RF(A,T,Δf)
        r2 = RF(A,T,Δf,0.0,r1.center)
        @test r1 ≈ r2
        @test isapprox(RF(A,T,Δf,0.0,r1.center,0.0), RF(A,T,Δf,0.0,r1.center,2π - 1e-5); rtol=1e-4)
        @test !isapprox(RF(A,T,Δf,0.0,r1.center,0.0), RF(A,T,Δf,0.0,r1.center,2π - 1e-5); rtol=1e-6)
        # Just checking to ensure that show() doesn't get stuck and that it is covered
        show(IOBuffer(), "text/plain", r1)
        @test true

        # Test Grad operations
        B1x, B1y, T = rand(3)
        A = B1x + im*B1y
        α = Complex(rand())
        rf = RF(A, T)
        rft = α * rf
        @test size(rf, 1) == 1
        @test rft.A ≈ α * rf.A
        rf_off = RF(0.0, T)
        rf_off_scaled = im * rf_off
        @test rf_off_scaled ≈ rf_off
        @test rf_off_scaled !== rf_off
        @test dur(rf) ≈ rf.T
        @test area(RF([0.0, 1.0, 0.0], [1.0, 1.0])) ≈ 1.0
        B1x, B1y, B2x, B2y, B3x, B3y, T1, T2, T3 = rand(9)
        rf1, rf2, rf3 = RF(B1x + im*B1y, T1), RF(B1x + im*B1y, T2), RF(B3x + im*B3y, T3)
        rv = [rf1; rf2; rf3 ;;]
        @test dur(rv) ≈ maximum(dur.(rv); dims=1)

        # Get Char from RF use and vice versa
        uses = [Excitation(), Refocusing(), Inversion(), Saturation(), Preparation(), Other(), Undefined()]
        chars = ['e', 'r', 'i', 's', 'p', 'o', 'u']
        for (use, char) in zip(uses, chars)
            @test KomaMRIBase.get_char_from_RF_use(use) == char
            @test KomaMRIBase.get_RF_use_from_char(Val(char)) == use
        end

        # Constant frequency offset: ψ is centered so the RF center has zero frame phase.
        rf = RF([1.0, 2.0, 1.0] .* 1e-6, 1e-3, 1000.0, 0.0; ϕ=0.3)
        @test KomaMRIBase.rf_frame_phase(rf) ≈ [0, π, -π, 0]

        # Frequency-modulated RF: ψ is generated from the integrated Δf waveform and centered.
        rf_fm = RF([1.0, 2.0, 1.0] .* 1e-6, 1e-3, [0.0, 1000.0, 0.0], 0.0)
        ψ = KomaMRIBase.rf_frame_phase(rf_fm)
        @test ψ[argmin(abs.(freq_times(rf_fm) .- (rf_fm.delay + rf_fm.center)))] ≈ 0

        # Discretization carries ψ so simulators can handle split RF blocks independently.
        seq = Sequence()
        seq += rf
        seqd = discretize(seq; sampling_params=Dict{String,Any}("Δt"=>1e-4, "Δt_rf"=>1e-4))
        center_idx = findfirst(==(rf.delay + rf.center), seqd.t)
        @test hasproperty(seqd, :ψ)
        @test !isnothing(center_idx)
        @test iszero(seqd.ψ[center_idx])
        @test length(seqd.ψ) == length(seqd.t)
        @test !all(iszero, seqd.ψ)
    end

    @testset "Delay" begin

        # Test delay construction
        T = 1e-3
        delay = Delay(T)
        duration = Duration(T)
        @test delay.T ≈ T
        @test dur(delay) ≈ T
        @test duration.T ≈ T
        @test dur(duration) ≈ T

        # Test delay construction error for negative values
        err = Nothing
        try Delay(-T) catch err end
        @test err isa ErrorException
        err = Nothing
        try Duration(-T) catch err end
        @test err isa ErrorException

        # Just checking to ensure that show() doesn't get stuck and that it is covered
        show(IOBuffer(), "text/plain", delay)
        show(IOBuffer(), "text/plain", duration)
        @test true

        # Test delay blocks
        seq = Sequence([Grad(0.0, 0.0)])
        ds = delay + seq
        @test dur(ds[1]) ≈ delay.T && dur(ds[2]) ≈ .0
        sd = seq + delay
        addblock!(sd, delay)
        @test dur(sd[1]) ≈ .0 && dur(sd[2]) ≈ delay.T && dur(sd[3]) ≈ delay.T
        sd = seq + duration
        addblock!(sd, duration)
        @test dur(sd[1]) ≈ .0 && dur(sd[2]) ≈ duration.T && dur(sd[3]) ≈ duration.T

        rf = RF(1e-6, 1e-3)
        block = Sequence()
        addblock!(block, rf, Delay(5e-3))
        @test length(block) == 1
        @test dur(block[1]) ≈ 5e-3

        block = Sequence()
        addblock!(block, rf, Delay(5e-4))
        @test length(block) == 1
        @test dur(block[1]) ≈ dur(rf)

        block = Sequence()
        addblock!(block, rf, Duration(5e-3))
        @test length(block) == 1
        @test dur(block[1]) ≈ 5e-3

        @test_throws ErrorException addblock!(Sequence(), rf, Duration(5e-4))
        @test_throws ErrorException addblock!(Sequence(), rf, Duration(5e-3), Duration(6e-3))

        adc = ADC(4, 2e-3)
        block = Sequence()
        addblock!(block, Delay(5e-4), rf, adc)
        @test length(block) == 1
        @test dur(block[1]) ≈ dur(adc)

        trigger = Trigger(1, 1, 1e-3, 2e-3)
        block = Sequence()
        addblock!(block, rf, trigger)
        @test length(block) == 1
        @test dur(block[1]) ≈ dur(trigger)

        block = Sequence() + (rf, Delay(5e-3))
        @test length(block) == 1
        @test dur(block[1]) ≈ 5e-3
        block = Sequence() + (rf, Duration(5e-3))
        @test length(block) == 1
        @test dur(block[1]) ≈ 5e-3

        block = Sequence()
        @addblock block += (rf, LabelSet(1, "LIN"), Duration(5e-3))
        @test length(block) == 1
        @test dur(block[1]) ≈ 5e-3
        @test block.EXT[1] == [LabelSet(1, "LIN")]

        t_delay, Δt_delay = KomaMRIBase.get_variable_times(Sequence() + delay)
        @test t_delay[2] ≈ 0.0
        @test t_delay[end-1] ≈ delay.T
        @test sum(Δt_delay) ≈ delay.T + 2KomaMRIBase.MIN_RISE_TIME

    end
    @testset "ADC" begin

        # Test ADC construction
        N, T, delay, Δf, ϕ  = 64, 1e-3, 2e-3, 1e-6, .25*π
        adc = ADC(N, T, delay, Δf, ϕ)

        adc1, adc2 = ADC(N, T), ADC(N,T,0,0,0)
        @test adc1 ≈ adc2

        adc1, adc2 = ADC(N, T, delay), ADC(N, T, delay, 0, 0)
        @test adc1 ≈ adc2

        adc1, adc2 = ADC(N, T, delay, Δf, ϕ), ADC(N, T, delay, Δf, ϕ)
        @test adc1 ≈ adc2
        @test isapprox(ADC(N, T, delay, Δf, 0.0), ADC(N, T, delay, Δf, 2π - 1e-5); rtol=1e-4)
        @test !isapprox(ADC(N, T, delay, Δf, 0.0), ADC(N, T, delay, Δf, 2π - 1e-5); rtol=1e-6)
        @test im * adc ≈ ADC(N, T, delay, Δf, ϕ + π / 2)
        @test adc * complex(-2.0) ≈ ADC(N, T, delay, Δf, mod(ϕ + π, 2π))
        adc_off = ADC(0, T)
        adc_off_scaled = im * adc_off
        @test adc_off_scaled ≈ adc_off
        @test adc_off_scaled !== adc_off

        # Test ADC construction errors for negative values
        err = Nothing
        try ADC(N, -T) catch err end
        @test err isa ErrorException
        try ADC(N, -T,  delay) catch err end
        @test err isa ErrorException
        try ADC(N,  T, -delay) catch err end
        @test err isa ErrorException
        try ADC(N, -T, -delay) catch err end
        @test err isa ErrorException
        try ADC(N, -T,  delay, Δf, ϕ) catch err end
        @test err isa ErrorException
        try ADC(N,  T, -delay, Δf, ϕ) catch err end
        @test err isa ErrorException
        try ADC(N, -T, -delay, Δf, ϕ) catch err end
        @test err isa ErrorException

        # Test ADC getproperties
        Nb, Tb, delayb, Δfb, ϕb  = 128, 2e-3, 4e-3, 2e-6, .125*π
        adb = ADC(Nb, Tb, delayb, Δfb, ϕb)
        adcs = [adc, adb]
        @test adcs.N ≈ [adc.N, adb.N] && adcs.T ≈ [adc.T, adb.T] && adcs.delay ≈ [adc.delay, adb.delay]
        @test adcs.Δf ≈ [adc.Δf, adb.Δf] && adcs.ϕ ≈ [adc.ϕ, adb.ϕ] && adcs.dur ≈ [adc.T + adc.delay, adb.T + adb.delay]

    end

    @testset "EXT" begin

        lInc = LabelInc(1,"LIN")
        lSet = LabelSet(1,"ECO")
        lSet2 = LabelSet(0,"LIN")
        lSetTRID = LabelSet(4,"TRID")
        trig = Trigger(0,1,100e-6,500e-6)
        @test dur(lInc) == 0.0
        @test dur(lSet) == 0.0
        @test dur(trig) == trig.delay + trig.duration

        d = Sequence([Grad(0,0.1)])
        seq = Sequence()
        d.EXT = [[lInc]]; 
        append!(seq, d)
        append!(seq, d)
        d.EXT = [[lInc,lSet]]
        append!(seq, d)
        d.EXT = [[lInc]]; 
        append!(seq, d)
        d.EXT = [[lSet2,trig]]; 
        append!(seq, d)
        d.EXT = [[]]; 
        append!(seq, d)

        @test seq.EXT[5][2] == trig && seq.EXT[5][1] == lSet2

        l = get_labels(seq)
        LIN_vec = [l[i].LIN for i in eachindex(l)] 
        @test LIN_vec == vec([1 2 3 4 0 0])

        ECO_vec = [l[i].ECO for i in eachindex(l)] 
        @test ECO_vec == vec([0 0 1 1 1 1])

        # Modification of the label directly in the sequence
        lSetPhs = LabelSet(2,"PHS")
        seq.EXT[4] = [lSetPhs]
        l = get_labels(seq)

        LIN_vec = [l[i].LIN for i in eachindex(l)] 
        @test LIN_vec == vec([1 2 3 3 0 0])
        PHS_vec = [l[i].PHS for i in eachindex(l)] 
        @test PHS_vec == vec([0 0 0 2 2 2])

        seq.EXT[6] = [lSetTRID]
        l = get_labels(seq)
        TRID_vec = [l[i].TRID for i in eachindex(l)]
        @test TRID_vec == vec([0 0 0 0 0 4])

    end

    @testset "DiscreteSequence" begin
        seq = PulseDesigner.EPI_example()
        sampling_params = KomaMRIBase.default_sampling_params()
        t, Δt = KomaMRIBase.get_variable_times(seq; Δt=sampling_params["Δt"], Δt_rf=sampling_params["Δt_rf"])
        seqd = KomaMRIBase.discretize(seq)
        i1, i2 = rand(1:Int(floor(0.5*length(seqd)))), rand(Int(ceil(0.5*length(seqd))):length(seqd))
        @test seqd[i1].t ≈ [t[i1]]
        @test seqd[i1:i2].t ≈ t[i1:i2]
        @test seqd[i1:i2].Δt ≈ Δt[i1:i2-1]

        T, N = 1.0, 4
        seq = Sequence()
        @addblock begin
            seq += RF(1.0e-6, 1.0)
            seq += Sequence([Grad(1.0e-3, 1.0)])
            seq += ADC(N, 1.0)
        end
        sampling_params = KomaMRIBase.default_sampling_params()
        sampling_params["Δt"], sampling_params["Δt_rf"] = T/N, T/N
        seqd1 = KomaMRIBase.discretize(seq[1]; sampling_params)
        seqd2 = KomaMRIBase.discretize(seq[2]; sampling_params)
        seqd3 = KomaMRIBase.discretize(seq[3]; sampling_params)
        # Block 1
        @test is_RF_on(seq[1]) == is_RF_on(seqd1)
        @test is_GR_on(seq[1]) == is_GR_on(seqd1)
        @test is_ADC_on(seq[1]) == is_ADC_on(seqd1)
        # Block 2
        @test is_RF_on(seq[2]) == is_RF_on(seqd2)
        @test is_GR_on(seq[2]) == is_GR_on(seqd2)
        @test is_ADC_on(seq[2]) == is_ADC_on(seqd2)
        # Block 3
        @test is_RF_on(seq[3]) == is_RF_on(seqd3)
        @test is_GR_on(seq[3]) == is_GR_on(seqd3)
        @test is_ADC_on(seq[3]) == is_ADC_on(seqd3)
        @test KomaMRIBase.is_GR_off(seqd) ==  !KomaMRIBase.is_GR_on(seqd)
        @test KomaMRIBase.is_RF_off(seqd) ==  !KomaMRIBase.is_RF_on(seqd)
        @test KomaMRIBase.is_ADC_off(seqd) == !KomaMRIBase.is_ADC_on(seqd)
    end

    @testset "large-time MRI event time-step collapse" begin
        T, offset = 1e-3, 200.0
        B1, Gx = 10e-6, 1e-3
        seq = Sequence()
        seq += Delay(offset)
        @addblock seq += RF(B1, T)
        seqd = KomaMRIBase.discretize(seq)
        area = KomaMRIBase.trapz(seqd.Δt, real.(seqd.B1))
        @test area ≈ B1 * T

        seq = Sequence()
        seq += Delay(offset)
        @addblock seq += Grad(Gx, T)
        seqd = KomaMRIBase.discretize(seq)
        area = KomaMRIBase.trapz(seqd.Δt, seqd.Gx)
        @test area ≈ Gx * T

        seq = Sequence()
        seq += Delay(offset)
        @addblock seq += ADC(2, T)
        seqd = KomaMRIBase.discretize(seq)
        @test all(diff(seqd.t) .> 0)
    end

    @testset "closing knot survives FP rebasing at large t0" begin
        for t0 in (270.0, 1500.0), event in (RF(1e-5, 1e-3), Grad(1e-3, 1e-3))
            ts = KomaMRIBase._reseparate_closing_knot!(t0 .+ KomaMRIBase.times(event))
            @test ts[end-1] < ts[end]
        end
    end

    @testset "get_samples and get_variable_times keep knots separated at large t0" begin
        T, offset = 1e-3, 300.0
        rf_seq = Sequence(); rf_seq += Delay(offset); @addblock rf_seq += RF(1e-5, T)
        gr_seq = Sequence(); gr_seq += Delay(offset); @addblock gr_seq += Grad(1e-3, T)
        ref_seq = Sequence(); ref_seq += Delay(0.5); @addblock ref_seq += Grad(1e-3, T)
        rf_t = KomaMRIBase.get_samples(rf_seq).rf.t
        @test rf_t[end-1] < rf_t[end]
        @test all(diff(KomaMRIBase.get_samples(gr_seq).gx.t) .> 0)
        @test length(KomaMRIBase.get_variable_times(gr_seq)[1]) ==
              length(KomaMRIBase.get_variable_times(ref_seq)[1])
    end

     @testset "SequenceFunctions" begin
        seq = PulseDesigner.EPI_example()
        t, Δt = KomaMRIBase.get_variable_times(seq; Δt=1)
        t_adc =  KomaMRIBase.get_adc_sampling_times(seq)
        M2, M2_adc = KomaMRIBase.get_slew_rate(seq)
        M2eddy, M2eddy_adc = KomaMRIBase.get_eddy_currents(seq)
        Gx, Gy, Gz = KomaMRIBase.get_grads(seq, t)
        Gmx, Gmy, Gmz = KomaMRIBase.get_grads(seq, reshape(t, 1, :))
        @test reshape(Gmx, :, 1) ≈ Gx && reshape(Gmy, :, 1) ≈ Gy && reshape(Gmz, :, 1) ≈ Gz
        @test is_ADC_on(seq) == is_ADC_on(seq, t)
        @test is_RF_on(seq) == is_RF_on(seq, t)
        @test KomaMRIBase.is_Delay(seq) == !(is_GR_on(seq) || is_RF_on(seq) || is_ADC_on(seq))
        @test size(M2, 1) == length(Δt) && size(M2_adc, 1) == length(t_adc)
        @test size(M2eddy, 1) == length(Δt) && size(M2eddy_adc, 1) == length(t_adc)

        # Just checking to ensure that show() doesn't get stuck and that it is covered
        show(IOBuffer(), "text/plain", seq)
        @test true

        α = rand()
        c = α + im*rand()
        x = copy(seq)
        x0 = copy(x)
        y = PulseDesigner.EPI_example()
        z = x + y
        @test z !== x
        @test z.GR.A ≈ [x0.GR  y.GR].A && z.RF.A ≈ [x0.RF y.RF].A && z.ADC.N ≈ [x0.ADC; y.ADC].N
        @test x ≈ x0
        z = x0 - y
        @test z.GR.A ≈ [x0.GR -y.GR].A
        z = -x0
        @test z.GR.A ≈ -x0.GR.A
        x_real = copy(x0)
        z = α * x_real
        @test z !== x_real
        @test z.GR.A ≈ α*x0.GR.A
        @test x_real.GR.A ≈ x0.GR.A
        rf_phase_after(c, rf) = is_RF_on(rf) ? mod(rf.ϕ + angle(c), 2π) : rf.ϕ
        adc_phase_after(c, adc) = is_ADC_on(adc) ? mod(adc.ϕ + angle(c), 2π) : adc.ϕ
        x_complex = copy(x0)
        z = c * x_complex
        @test z !== x_complex
        @test all(rfz.A ≈ abs(c) * rfx.A &&
                  rfz.ϕ ≈ rf_phase_after(c, rfx)
                  for (rfz, rfx) in zip(z.RF, x0.RF))
        @test all(adcz.ϕ ≈ adc_phase_after(c, adcx) for (adcz, adcx) in zip(z.ADC, x0.ADC))
        @test x_complex.GR.A ≈ x0.GR.A
        @test_throws MethodError x0 * α
        @test_throws MethodError x0 * c
        @test_throws MethodError x0 * rotz(π / 4)
        z = x0 / α
        @test z.GR.A ≈ x0.GR.A/α
        @test size(y) == size(y.GR[1,:])
        z = x0 + x0.GR[3,1]
        @test z.GR.A[1, end] ≈ x0.GR[3,1].A
        z = x0.GR[3,1] + x0
        @test z.GR.A[1, 1] ≈ x0.GR[3,1].A
        z = x + x.RF[1,1]
        @test z.RF.A[1, end] ≈ x.RF[1,1].A
        z = x.RF[1,1] + x
        @test z.RF.A[1, 1] ≈ x.RF[1,1].A
        z = x + x.ADC[3,1]
        @test z.ADC.N[end] ≈ x.ADC[3,1].N
        z = x.ADC[3,1] + x
        @test z.ADC.N[1] ≈ x.ADC[3,1].N

        seq_alias = Sequence(
            reshape([Grad(0.01, 0.01, 0.01), Grad(0.02, 0.01, 0.01)], 1, :),
            reshape([RF(1e-6, 1e-3), RF(2e-6, 1e-3)], 1, :),
            [ADC(4, 1e-3), ADC(8, 2e-3)],
            [0.02, 0.03],
            [Extension[], Extension[]],
            Dict{String, Any}(),
        )
        seq_alias[2].GR[1,1].A = 1e-3
        @test seq_alias.GR[1, 2].A == 1e-3
        seq_alias[2].RF[1].A = 3e-6
        @test seq_alias.RF[1, 2].A == 3e-6 + 0im
        seq_alias[2].ADC[1].N = 16
        @test seq_alias.ADC[2].N == 16
        seq_alias[2].DUR[1] = 0.05
        @test seq_alias.DUR[2] == 0.05
        seq_alias[1].GR[1,1] = Grad(0.5, 0.01, 0.01)
        @test seq_alias.GR[1, 1].A == 0.5
        seq_alias[1:2].DUR .= [0.06, 0.07]
        @test seq_alias.DUR[1] == 0.06
        @test seq_alias.DUR[2] == 0.07
        seq_alias[BitVector([false, true])].ADC[1].N = 32
        @test seq_alias.ADC[2].N == 32
    end
    @testset "Sequence comparison" begin
        seq1 = PulseDesigner.EPI_example()
        seq2 = deepcopy(seq1)
        @test seq1 ≈ seq2 # Basic equality check
        rf = seq1[1].RF[1]
        seq2 = Sequence([
            Sequence(
            seq1.GR[:, 1],
            [RF([rf.A, rf.A, rf.A], rf.T, rf.Δf, rf.delay, rf.center, rf.use)],
            seq1.ADC[1],
            seq1.DUR[1],
            seq1.EXT[1],
            seq1.DEF,
            ),
            seq1[2:end],
        ])
        @test seq1 ≉ seq2 # Structural equality does not treat different RF samplings as equal
        seq2.EXT[1] = [LabelInc(2, "LIN")]
        @test seq1 ≉ seq2 # Check that sequences are not equal because one of them has extensions
        seq1.EXT[1] = [LabelSet(2, "ECO")]
        @test seq1 ≉ seq2 # Check that sequences are not equal because they have different extensions
        seq2 = deepcopy(seq1)
        seq1.EXT[1] = [LabelInc(2, "LIN"), LabelSet(2, "ECO")]
        seq2.EXT[1] = [LabelInc(2, "LIN"), LabelSet(2, "ECO")]
        @test seq1 ≈ seq2 # Check that sequences are equal because they have the same extensions
    end
    @testset "Quaternion rotations" begin
        @testset "Quaternion conversion" begin
            q = QuaternionRot(0, 0, 0, 1)
            @test KomaMRIBase.rotation_matrix(q) == [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]
            @test KomaMRIBase.rotation_matrix(QuaternionRot([1 0 0; 0 0 -1; 0 1 0])) ≈ [1 0 0; 0 0 -1; 0 1 0]
            @test KomaMRIBase.rotation_matrix(QuaternionRot(rotz(π / 3))) ≈ rotz(π / 3)
            @test_throws DimensionMismatch QuaternionRot([1 0; 0 1])
            @test_throws ArgumentError QuaternionRot(zeros(3, 3))
        end

        @testset "Apply rotations keeps extension metadata" begin
            seq = Sequence()
            @addblock seq += (QuaternionRot(0, 0, 0, 1), LabelSet(1, "LIN"), x=Grad(1e-3, 1e-3), y=Grad(2e-3, 1e-3))
            seq.DEF["RequiredExtensions"] = ["ROTATIONS", "LABELSET"]
            rotated = apply_rotations(seq)
            @test rotated !== seq
            @test rotated.GR[1, 1].A == -1e-3
            @test rotated.GR[2, 1].A == -2e-3
            @test rotated.GR[3, 1].A == 0.0
            @test rotated.EXT[1] == [QuaternionRot(0, 0, 0, 1), LabelSet(1, "LIN")]
            @test rotated.DEF["RequiredExtensions"] == ["ROTATIONS", "LABELSET"]
            @test seq.GR[1, 1].A == 1e-3
            @test any(ext -> ext isa QuaternionRot, seq.EXT[1])

            mixed = Sequence([Grad(1e-3, 1e-3); Grad(2e-3, 1e-3, 0.0, 0.0, 0.2e-3); Grad(0.0, 0.0);;])
            mixed.EXT[1] = [QuaternionRot(cos(π / 8), 0, 0, sin(π / 8))]
            mixed_rotated = apply_rotations(mixed)
            @test mixed_rotated !== mixed
            @test any(ext -> ext isa QuaternionRot, mixed_rotated.EXT[1])
        end
    end
    @testset "Check Scanner Constraints" begin
        sys = Scanner()
        seq = PulseDesigner.EPI_example(; sys)
        @test isnothing(check_hw_limits(seq, sys))
    end
    @testset "Sequence timing and hardware metadata" begin
        seq = Sequence()
        @test seq.DEF["BlockDurationRaster"] == 10e-6
        @test seq.DEF["GradientRasterTime"] == 10e-6
        @test seq.DEF["RadiofrequencyRasterTime"] == 1e-6
        @test seq.DEF["AdcRasterTime"] == 100e-9
        @test isinf(seq.DEF["MaxB1"])
        @test isinf(seq.DEF["MaxGrad"])
        @test isinf(seq.DEF["MaxSlew"])

        sys = Scanner(B1=20e-6, Gmax=40e-3, Smax=150.0, ADC_Δt=2e-6, DUR_Δt=20e-6, GR_Δt=10e-6, RF_Δt=2e-6)
        seq = Sequence(sys)
        @test seq.DEF["BlockDurationRaster"] == sys.DUR_Δt
        @test seq.DEF["AdcRasterTime"] == sys.ADC_Δt
        @test seq.DEF["MaxB1"] == sys.B1
        @test seq.DEF["MaxGrad"] == sys.Gmax

        file_def = KomaMRIBase._sequence_def_from_pulseq(Dict("BlockDurationRaster" => 2e-6))
        @test file_def["BlockDurationRaster"] == 2e-6
        @test isinf(file_def["MaxGrad"])

        @testset "Pulseq checkTiming raster checks" begin
            block_raster = 10e-6
            rf_raster = 1e-6
            adc_raster = 100e-9
            seq = Sequence()
            @addblock seq += (RF(1e-6, block_raster), x=Grad(1e-3, block_raster))
            @test isnothing(check_timing(seq))
            @test isnothing(check_hw_limits(seq))
            @test_throws ErrorException check_hw_limits(seq, Scanner(B1=0.5e-6))

            compact_rf = Sequence()
            @addblock compact_rf += (RF([1e-6, 1e-6], rf_raster, 0.0, rf_raster / 2), Duration(block_raster))
            @test isnothing(check_timing(compact_rf))

            adc_samples = 2
            adc_dwell = adc_raster
            compact_adc = Sequence()
            @addblock compact_adc += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_dwell / 2), Duration(block_raster))
            @test isnothing(check_timing(compact_adc))

            adc_delay_raster = block_raster
            raster_sys = Scanner(RF_Δt=adc_delay_raster, ADC_Δt=adc_raster, ADC_dead_time=0.0)
            adc_raster_seq = Sequence(raster_sys)
            adc_samples = 111
            adc_dwell = 44 * adc_raster
            adc_pulseq_delay = 15 * adc_delay_raster
            adc_delay = adc_pulseq_delay + adc_dwell / 2
            adc_duration = ceil((adc_pulseq_delay + adc_samples * adc_dwell) / raster_sys.DUR_Δt) * raster_sys.DUR_Δt
            @addblock adc_raster_seq += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_delay), Duration(adc_duration))
            @test isnothing(check_timing(adc_raster_seq))

            adc_bad_dwell = Sequence(Scanner(ADC_Δt=adc_raster))
            bad_dwell = adc_raster + adc_raster / 2
            @addblock adc_bad_dwell += (ADC(2, bad_dwell, bad_dwell / 2), Duration(block_raster))
            @test_throws ErrorException check_timing(adc_bad_dwell)

            adc_bad_delay = Sequence(Scanner(RF_Δt=adc_delay_raster, ADC_Δt=adc_raster))
            bad_pulseq_delay = adc_delay_raster + adc_raster
            @addblock adc_bad_delay += (ADC(2, adc_raster, bad_pulseq_delay + adc_raster / 2), Duration(2block_raster))
            @test_throws ErrorException check_timing(adc_bad_delay)

            seq_bad = Sequence()
            @addblock seq_bad += RF(1e-6, block_raster + rf_raster / 2)
            @test_throws ErrorException check_timing(seq_bad)
        end

        @testset "Pulseq shapes-and-times duration checks" begin
            block_raster = 10e-6
            rf_raster = 1e-6
            adc_raster = 100e-9
            event_duration = 2 * block_raster
            adc_edge = Sequence(Scanner(
                DUR_Δt=adc_raster, RF_Δt=adc_raster, ADC_Δt=adc_raster,
                ADC_dead_time=0.0,
            ))
            @addblock adc_edge += ADC(2, adc_raster, adc_raster / 2)
            @test dur(adc_edge) ≈ 3adc_raster / 2
            @test_throws ErrorException check_timing(adc_edge)
            adc_edge.DUR[1] = 2adc_raster
            @test isnothing(check_timing(adc_edge))

            rf_too_long = Sequence()
            @addblock rf_too_long += RF(1e-6, event_duration)
            rf_too_long.DUR[1] = block_raster
            @test_throws ErrorException check_timing(rf_too_long)

            grad_too_long = Sequence()
            @addblock grad_too_long += (x=Grad(1e-3, event_duration))
            grad_too_long.DUR[1] = block_raster
            @test_throws ErrorException check_timing(grad_too_long)

            adc_too_long = Sequence(Scanner(DUR_Δt=adc_raster, RF_Δt=adc_raster, ADC_Δt=adc_raster))
            @addblock adc_too_long += ADC(2, adc_raster, adc_raster / 2)
            adc_too_long.DUR[1] = adc_raster
            @test_throws ErrorException check_timing(adc_too_long)

            trigger_too_long = Sequence()
            @addblock trigger_too_long += Trigger(1, 1, 0.0, event_duration)
            trigger_too_long.DUR[1] = block_raster
            @test_throws ErrorException check_timing(trigger_too_long)
        end

        @testset "Pulseq checkTiming dead-time and ring-down checks" begin
            rf_raster = 1e-6
            adc_raster = 100e-9
            rf_dead_time = 10e-6
            rf_ring_down_time = 20e-6
            adc_dead_time = 10e-6
            deadtime_sys = Scanner(B1=Inf, Gmax=Inf, Smax=Inf, ADC_Δt=adc_raster, RF_Δt=rf_raster, RF_dead_time=rf_dead_time, RF_ring_down_time=rf_ring_down_time, ADC_dead_time=adc_dead_time)
            rf_delay = rf_dead_time + rf_raster / 2
            rf_duration = rf_raster
            rf_end = rf_dead_time + 2rf_raster
            rf_block_duration = ceil((rf_end + rf_ring_down_time) / deadtime_sys.DUR_Δt) * deadtime_sys.DUR_Δt
            rf_deadtime_ok = Sequence(deadtime_sys)
            @addblock rf_deadtime_ok += (RF([1e-6, 1e-6], rf_duration, 0.0, rf_delay), Duration(rf_block_duration))
            @test isnothing(check_timing(rf_deadtime_ok, deadtime_sys))

            rf_deadtime_bad = Sequence(deadtime_sys)
            @addblock rf_deadtime_bad += (RF([1e-6, 1e-6], rf_duration, 0.0, rf_dead_time), Duration(rf_block_duration))
            @test_throws ErrorException check_timing(rf_deadtime_bad, deadtime_sys)

            rf_ringdown_bad = Sequence(deadtime_sys)
            @addblock rf_ringdown_bad += (RF([1e-6, 1e-6], rf_duration, 0.0, rf_delay), Duration(rf_block_duration - deadtime_sys.DUR_Δt))
            @test_throws ErrorException check_timing(rf_ringdown_bad, deadtime_sys)

            adc_samples = 2
            adc_dwell = adc_raster
            adc_delay = adc_dead_time + adc_dwell / 2
            adc_end_with_dead_time = adc_dead_time + adc_samples * adc_dwell + adc_dead_time
            adc_block_duration = ceil(adc_end_with_dead_time / deadtime_sys.DUR_Δt) * deadtime_sys.DUR_Δt
            adc_deadtime_ok = Sequence(deadtime_sys)
            @addblock adc_deadtime_ok += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_delay), Duration(adc_block_duration))
            @test isnothing(check_timing(adc_deadtime_ok, deadtime_sys))

            adc_deadtime_bad = Sequence(deadtime_sys)
            @addblock adc_deadtime_bad += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_delay - adc_dead_time / 10), Duration(adc_block_duration))
            @test_throws ErrorException check_timing(adc_deadtime_bad, deadtime_sys)

            adc_post_deadtime_bad = Sequence(deadtime_sys)
            @addblock adc_post_deadtime_bad += (ADC(adc_samples, (adc_samples - 1) * adc_dwell, adc_delay), Duration(adc_block_duration - deadtime_sys.DUR_Δt))
            @test_throws ErrorException check_timing(adc_post_deadtime_bad, deadtime_sys)
        end

        @testset "C++ interpreter extension constraints" begin
            block_raster = 10e-6
            duplicate_trigger = Sequence()
            @addblock duplicate_trigger += (Trigger(1, 1, 0.0, 0.0), Trigger(1, 2, 0.0, 0.0), Duration(block_raster))
            @test_throws ErrorException check_timing(duplicate_trigger)

            duplicate_rotation = Sequence()
            @addblock duplicate_rotation += (QuaternionRot(1, 0, 0, 0), QuaternionRot(1, 0, 0, 0), Duration(block_raster))
            @test_throws ErrorException check_timing(duplicate_rotation)
        end
    end
end

@testitem "PulseDesigner" tags=[:base] begin
    import KomaMRIBase: rotation_matrix
    using Unitful

    @testset "RF_sinc" begin
        sys = Scanner()
        B1 = 23.4e-6 # For 90 deg flip angle
        Trf = 1e-3
        rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4)
        @test round(get_flip_angles(rf)[1]) ≈ 90
    end
    @testset "build_block_pulse" begin
        sys = Scanner(B0=3u"T")
        @test sys.B0 == 3.0

        sys = Scanner(B1=γ * 20e-6 * u"Hz")
        @test sys.B1 ≈ 20e-6

        sys = Scanner(B1=Inf, RF_Δt=1u"μs", RF_dead_time=0u"s", RF_ring_down_time=0u"s")
        seq = PulseDesigner.build_block_pulse(π * u"rad"; duration=1u"ms", sys)
        @test get_flip_angles(seq)[1] ≈ 180.0

        seq = PulseDesigner.build_block_pulse(90u"°"; duration=2u"ms", sys)
        @test get_flip_angles(seq)[1] ≈ 90.0

        sys = Scanner(B1=Inf, RF_Δt=1u"μs", DUR_Δt=10u"μs", RF_dead_time=100u"μs", RF_ring_down_time=30u"μs")
        seq = PulseDesigner.build_block_pulse(8u"°"; duration=200u"μs", sys)
        @test dur(seq) ≈ 330e-6

        rf = PulseDesigner.build_block_pulse(π * u"rad"; bandwidth=0.3u"kHz", sys).RF[1, 1]
        @test rf.T == round(1 / (4 * 0.3e3) / sys.RF_Δt) * sys.RF_Δt

        seq = PulseDesigner.build_block_pulse(π * u"rad"; bandwidth=1u"kHz", time_bw_product=5, sys)
        rf = seq.RF[1, 1]
        @test rf.T == 5e-3

        seq = PulseDesigner.build_block_pulse(
            90u"°"; duration=1u"ms", freq_offset=1u"kHz", phase_offset=90u"°",
            sys, use=Excitation(),
        )
        rf = seq.RF[1, 1]
        @test rf.Δf == 1e3
        @test rf.ϕ == π / 2

    end
    @testset "Unitful API parity" begin
        # Pulseq-default plain gradient numerics should match equivalent Unitful SI inputs.
        @test (2.0u"ms" |> to_SI) == 2e-3
        @test ([256.0u"mm", 5.0u"mm"] .|> to_SI) == [0.256, 0.005]
        @test PulseDesigner.ceil_to_raster(1.011u"ms", 10u"μs") == 1.02u"ms"
        @test PulseDesigner.floor_to_raster(1.019u"ms", 10u"μs") == 1.01u"ms"
        @test PulseDesigner.round_to_raster(1.014u"ms", 10u"μs") == 1.01u"ms"
        @test PulseDesigner.raster_samples(1.014u"ms", 10u"μs") == 101

        sys = Scanner(
            B1=Inf, Gmax=40e-3, Smax=170.0, ADC_Δt=2e-6, DUR_Δt=10e-6,
            GR_Δt=10e-6, RF_Δt=1e-6, RF_ring_down_time=0.0,
            RF_dead_time=0.0, ADC_dead_time=6e-6,
        )
        sysu = Scanner(
            B1=Inf, Gmax=40u"mT/m", Smax=170u"T/m/s", ADC_Δt=2u"μs",
            DUR_Δt=10u"μs", GR_Δt=10u"μs", RF_Δt=1u"μs",
            RF_ring_down_time=0u"s", RF_dead_time=0u"s", ADC_dead_time=6u"μs",
        )
        physical_area = 8e-6
        physical_amplitude = 1e-3
        @test PulseDesigner.build_trapezoid(:x; area=γ * physical_area, duration=1e-3, sys) ≈
            PulseDesigner.build_trapezoid(:x; area=8e-6u"T*s/m", duration=1u"ms", sys=sysu)
        @test PulseDesigner.build_arbitrary_grad(
            :x, γ .* [0, physical_amplitude, 0]; sys,
        ) ≈
            PulseDesigner.build_arbitrary_grad(:x, [0, 1, 0]u"mT/m"; sys=sysu)
        @test PulseDesigner.build_extended_trapezoid(
            :x, [0, 0.5e-3, 1e-3], γ .* [0, physical_amplitude, 0]; sys,
        ) ≈ PulseDesigner.build_extended_trapezoid(
            :x, [0, 0.5, 1]u"ms", [0, 1, 0]u"mT/m"; sys=sysu,
        )
        @test PulseDesigner.build_extended_trapezoid_area(
            :x, 0.0, 0.0, γ * 10e-6; sys,
        ) ≈
            PulseDesigner.build_extended_trapezoid_area(
                :x, 0u"mT/m", 0u"mT/m", 10e-6u"T*s/m"; sys=sysu,
            )
        @test PulseDesigner.build_block_pulse(π / 2; duration=1e-3, sys) ≈
            PulseDesigner.build_block_pulse(90u"°"; duration=1u"ms", sys=sysu)
        @test PulseDesigner.build_sinc_pulse(π / 2; duration=1e-3, sys) ≈
            PulseDesigner.build_sinc_pulse(90u"°"; duration=1u"ms", sys=sysu)
        @test PulseDesigner.build_arbitrary_rf([1, 2, 1], π / 2; dwell=2e-6, sys) ≈
            PulseDesigner.build_arbitrary_rf([1, 2, 1], 90u"°"; dwell=2u"μs", sys=sysu)
        @test PulseDesigner.build_gauss_pulse(π / 2; duration=1e-3, sys) ≈
            PulseDesigner.build_gauss_pulse(90u"°"; duration=1u"ms", sys=sysu)
        @test PulseDesigner.build_adiabatic_pulse(
            :wurst; duration=4e-3, dwell=2e-6, bandwidth=6e3, sys,
        ) ≈ PulseDesigner.build_adiabatic_pulse(
            :wurst; duration=4u"ms", dwell=2u"μs", bandwidth=6u"kHz", sys=sysu,
        )
        @test PulseDesigner.build_adc(4, 2e-6; sys) ≈
            PulseDesigner.build_adc(4, 2u"μs"; sys=sysu)
        @test PulseDesigner.build_delay(1e-3; sys) ≈
            PulseDesigner.build_delay(1u"ms"; sys=sysu)
        @test PulseDesigner.build_rotation(π / 2; sys) ≈
            PulseDesigner.build_rotation(90u"deg"; sys=sysu)
        @test PulseDesigner.build_trigger(:physio1; delay=20e-6, duration=100e-6, sys) ≈
            PulseDesigner.build_trigger(
                :physio1; delay=20u"μs", duration=100u"μs", sys=sysu,
            )
        @test PulseDesigner.build_digital_output_pulse(
            :osc0; delay=20e-6, duration=100e-6, sys,
        ) ≈ PulseDesigner.build_digital_output_pulse(
            :osc0; delay=20u"μs", duration=100u"μs", sys=sysu,
        )
    end
    @testset "build_trapezoid" begin
        # Cover the constructor branches that solve amplitude from area and preserve amplitude directly.
        sys = Scanner(Gmax=Inf, Smax=Inf, GR_Δt=10e-6)
        target_area = γ * 10e-6
        expected_area = target_area / γ
        duration = 1e-3
        rise_time = sys.GR_Δt
        grad = PulseDesigner.make_trapezoid(;
            area=target_area, duration, rise_time, sys,
        )
        @test area(grad) ≈ expected_area

        target_amplitude = γ * 20e-3
        expected_amplitude = target_amplitude / γ
        grad = PulseDesigner.make_trapezoid(;
            amplitude=target_amplitude, flat_time=0.8e-3, rise_time, sys,
        )
        @test grad.A ≈ expected_amplitude

        # Unitful gradient areas and Pulseq-style Hz/m amplitudes convert to SI T/m values.
        sys = Scanner(Gmax=40u"mT/m", Smax=170u"T/m/s", GR_Δt=10u"μs")
        seq = PulseDesigner.build_trapezoid(:x; area=11e-6u"T*s/m", sys)
        grad = only(seq.GR.x)
        @test area(grad) ≈ 1.1e-5

        seq = PulseDesigner.build_trapezoid(:x; amplitude=20u"mT/m", flat_time=0.8u"ms", sys)
        grad = only(seq.GR.x)
        @test grad.A ≈ 20e-3

        seq = PulseDesigner.build_trapezoid(:x; area=10e-6u"T*s/m", duration=1u"ms", sys)
        grad = only(seq.GR.x)
        @test area(grad) ≈ 1e-5

        seq = PulseDesigner.build_trapezoid(:x; flat_area=8e-6u"T*s/m", flat_time=0.8u"ms", sys)
        grad = only(seq.GR.x)
        @test grad.A ≈ 10e-3

        amplitude = γ * 20e-3 * u"Hz/m"
        seq = PulseDesigner.build_trapezoid(:x; amplitude, flat_time=0.8u"ms", sys)
        grad = only(seq.GR.x)
        @test grad.A ≈ 20e-3
    end
    @testset "build_extended_trapezoid" begin
        sys = Scanner(Gmax=40u"mT/m", Smax=170u"T/m/s", GR_Δt=10u"μs")
        times = [0.0, 0.5, 1.0]u"ms"
        amplitudes = [0.0, 20.0, 0.0]u"mT/m"
        seq = PulseDesigner.build_extended_trapezoid(:x, times, amplitudes; sys)
        grad = only(seq.GR.x)
        @test area(grad) ≈ 1e-5
    end
    @testset "build_arbitrary_grad" begin
        # Arbitrary gradients infer edge samples, oversampling delay, and reject invalid waveforms.
        sys = Scanner(Gmax=40e-3, Smax=170.0, GR_Δt=10e-6)
        physical_waveform = [0, 1e-3, 0]
        waveform = γ .* physical_waveform
        grad = PulseDesigner.make_arbitrary_grad(waveform; sys)
        @test grad.first ≈ -maximum(physical_waveform) / 2
        @test grad.last ≈ -maximum(physical_waveform) / 2
        @test grad.T == 2sys.GR_Δt

        oversampled_physical_waveform = [0, 0.2e-3, 0.4e-3, 0.2e-3, 0]
        oversampled_waveform = γ .* oversampled_physical_waveform
        oversampled_edge = -γ * 0.2e-3
        seq = PulseDesigner.build_arbitrary_grad(
            :y, oversampled_waveform;
            oversampling=true, first=oversampled_edge, last=oversampled_edge, sys,
        )
        grad = only(seq.GR.y)
        @test grad.T == 2sys.GR_Δt
        @test grad.rise == sys.GR_Δt / 2

        single_sample = γ * 0.5e-3
        grad = PulseDesigner.make_arbitrary_grad([single_sample]; first=0, last=0, sys)
        @test grad.T == 0.0
        @test grad.first == 0.0
        @test grad.last == 0.0

        @test_throws ErrorException PulseDesigner.make_arbitrary_grad([0.0]; sys)
        @test_throws ErrorException PulseDesigner.make_arbitrary_grad(
            γ .* [0, 1e-3]; oversampling=true, sys,
        )
        @test_throws ErrorException PulseDesigner.make_arbitrary_grad(
            γ .* [0, 50e-3, 0]; sys,
        )
        @test_throws ErrorException PulseDesigner.build_arbitrary_grad(
            :q, γ .* [0, 1e-3, 0]; sys,
        )
    end
    @testset "build_extended_trapezoid_area" begin
        sys = Scanner(Gmax=40u"mT/m", Smax=170u"T/m/s", GR_Δt=10u"μs")
        cases = (
            (0u"mT/m", 0u"mT/m", 0u"T*s/m"),
            (10u"mT/m", 10u"mT/m", 0u"T*s/m"),
            (0u"mT/m", 0u"mT/m", 10e-6u"T*s/m"),
            (10u"mT/m", 5u"mT/m", 2e-6u"T*s/m"),
            (0u"mT/m", 0u"mT/m", -5e-6u"T*s/m"),
        )
        for (grad_start, grad_end, target_area) in cases
            seq = PulseDesigner.build_extended_trapezoid_area(
                :x, grad_start, grad_end, target_area; sys,
            )
            grad = only(seq.GR.x)
            @test area(grad) ≈ ustrip(u"T*s/m", target_area)
            @test grad.first ≈ ustrip(u"T/m", grad_start)
            @test grad.last ≈ ustrip(u"T/m", grad_end)
        end
    end
    @testset "build_rotation" begin
        @test rotation_matrix(PulseDesigner.make_rotation(π / 2)) ≈ rotz(π / 2)
        @test rotation_matrix(PulseDesigner.make_rotation(π / 2, π / 4)) ≈
            rotz(π / 2) * roty(π / 4)
        @test rotation_matrix(PulseDesigner.make_rotation([1, 0, 0], π / 2)) ≈
            rotx(π / 2)
        @test rotation_matrix(PulseDesigner.make_rotation([0, -1, 0], π / 3)) ≈
            roty(-π / 3)
        @test rotation_matrix(PulseDesigner.make_rotation([0, 0, 0, 2])) ≈
            rotz(π)
        @test rotation_matrix(PulseDesigner.make_rotation(rotz(π / 3))) ≈
            rotz(π / 3)

        seq = PulseDesigner.build_rotation(90u"°")
        @test rotation_matrix(only(seq.EXT[1])) ≈ rotz(π / 2)
        @test rotation_matrix(PulseDesigner.make_rotation(90u"°", 30u"°")) ≈
            rotz(π / 2) * roty(π / 6)
        @test rotation_matrix(PulseDesigner.make_rotation(-90u"°")) ≈
            rotz(-π / 2)
        @test rotation_matrix(PulseDesigner.make_rotation([0, 1, 0], 45u"°")) ≈
            roty(π / 4)
        @test rotation_matrix(PulseDesigner.make_rotation([0, -1, 0], 45u"°")) ≈
            roty(-π / 4)
        @test_throws DimensionMismatch PulseDesigner.make_rotation(1u"ppm")
    end
    @testset "build_trigger" begin
        sys = Scanner(GR_Δt=10u"μs")
        seq = PulseDesigner.build_trigger(:physio2; delay=20u"μs", duration=100u"μs", sys)
        trigger = only(seq.EXT[1])
        @test trigger == Trigger(2, 2, 20e-6, 100e-6)
    end
    @testset "build_digital_output_pulse" begin
        sys = Scanner(GR_Δt=10u"μs")
        seq = PulseDesigner.build_digital_output_pulse(
            :ext1; delay=500u"μs", duration=100u"μs", sys,
        )
        output = only(seq.EXT[1])
        @test output == Trigger(1, 3, 500e-6, 100e-6)
    end
    @testset "build_label" begin
        seq = PulseDesigner.build_label(:SET, :LIN, true)
        append!(seq, PulseDesigner.build_label(:INC, :LIN, 2))

        labels = get_labels(seq)
        @test [label.LIN for label in labels] == [1, 3]
    end
    @testset "build_delay" begin
        seq = PulseDesigner.build_delay(2u"ms")
        @test only(seq.DUR) == 2e-3
    end
    @testset "Pulseq RF/ADC timing" begin
        sys = Scanner(
            RF_Δt=1e-6, RF_ring_down_time=30e-6,
            ADC_Δt=2e-6, ADC_dead_time=6e-6, DUR_Δt=10e-6,
        )

        # Uniform RF stores delay at the first sample; sys-aware timing writes the
        # Pulseq delay at the preceding sample edge and includes RF ringdown.
        rf_waveform = ComplexF64[1, 1, 1]
        rf_pulseq_delay = 3sys.RF_Δt
        rf_sample_offset = sys.RF_Δt / 2
        rf_duration = (length(rf_waveform) - 1) * sys.RF_Δt
        rf = RF(rf_waveform, rf_duration, 0.0, rf_pulseq_delay + rf_sample_offset; center=sys.RF_Δt)
        @test dur(rf) ≈ delay(rf) + rf_duration
        @test delay(rf) ≈ rf_pulseq_delay + rf_sample_offset
        @test rf_center(rf) ≈ sys.RF_Δt
        @test delay(rf, sys) ≈ rf_pulseq_delay
        @test rf_center(rf, sys) ≈ rf_center(rf) + rf_sample_offset
        @test dur(rf, sys) ≈ delay(rf, sys) + length(rf.A) * sys.RF_Δt + sys.RF_ring_down_time

        # ADC stores delay at the first sample; sys-aware timing writes the Pulseq
        # delay at dwell/2 before that sample and includes post-ADC dead time.
        adc_samples = 4
        adc_dwell = sys.ADC_Δt
        adc_pulseq_delay = 3sys.ADC_Δt
        adc = PulseDesigner.make_adc(adc_samples, adc_dwell; delay=adc_pulseq_delay - sys.ADC_dead_time, sys)
        @test dur(adc) ≈ delay(adc) + (adc_samples - 1) * adc_dwell
        @test delay(adc) ≈ adc_pulseq_delay + adc_dwell / 2
        @test delay(adc, sys) ≈ adc_pulseq_delay
        @test dur(adc, sys) ≈ adc_pulseq_delay + adc_samples * adc_dwell + sys.ADC_dead_time
        adc_block = PulseDesigner.build_adc(adc_samples, adc_dwell; delay=adc_pulseq_delay - sys.ADC_dead_time, sys)
        @test dur(adc_block, sys) ≈ dur(adc, sys)
        @test only(adc_block.DUR) ≈ PulseDesigner.ceil_to_raster(dur(adc, sys), sys.DUR_Δt)

        # Block-pulse RF has no sample-edge shift; sys-aware duration only adds ringdown.
        rf_block = PulseDesigner.build_block_pulse(90u"°"; duration=2u"μs", sys)
        @test dur(rf_block[1], sys) ≈ only(rf_block.DUR)
        block_rf_delay = sys.RF_dead_time
        block_rf_duration = 2sys.RF_Δt
        block_rf = RF(1e-6, block_rf_duration, 0.0, block_rf_delay; center=block_rf_duration / 2, use=Excitation())
        @test delay(block_rf, sys) ≈ delay(block_rf)
        @test rf_center(block_rf, sys) ≈ rf_center(block_rf)
        @test dur(block_rf, sys) ≈ delay(block_rf) + sum(block_rf.T) + sys.RF_ring_down_time

        # Compact RF timing uses type dispatch: default-raster RF writes id 0,
        # half-raster RF writes an explicit time shape with RF-raster duration.
        rf_samples = ComplexF64[1, 2, 1]
        uniform_rf = RF(rf_samples, (length(rf_samples) - 1) * sys.RF_Δt, 0.0, rf_pulseq_delay + sys.RF_Δt / 2; center=sys.RF_Δt, use=Excitation())
        @test dwell(uniform_rf, sys) ≈ sys.RF_Δt
        @test delay(uniform_rf) ≈ rf_pulseq_delay + dwell(uniform_rf) / 2
        @test delay(uniform_rf, sys) ≈ rf_pulseq_delay
        @test rf_center(uniform_rf) ≈ dwell(uniform_rf)
        @test rf_center(uniform_rf, sys) ≈ rf_center(uniform_rf) + dwell(uniform_rf) / 2
        @test dur(uniform_rf) ≈ delay(uniform_rf) + 2dwell(uniform_rf)
        @test dur(uniform_rf, sys) ≈ delay(uniform_rf, sys) + length(uniform_rf.A) * dwell(uniform_rf) + sys.RF_ring_down_time

        half_raster_dwell = sys.RF_Δt / 2
        half_raster_rf = RF(rf_samples, (length(rf_samples) - 1) * half_raster_dwell, 0.0, rf_pulseq_delay + half_raster_dwell / 2; center=half_raster_dwell, use=Excitation())
        @test dwell(half_raster_rf, sys) ≈ sys.RF_Δt / 2
        @test delay(half_raster_rf, sys) ≈ rf_pulseq_delay
        @test rf_center(half_raster_rf, sys) ≈ rf_center(half_raster_rf) + dwell(half_raster_rf) / 2
        half_raster_shape_duration = PulseDesigner.ceil_to_raster(
            dur(half_raster_rf) + dwell(half_raster_rf) / 2 - delay(half_raster_rf, sys),
            sys.RF_Δt,
        )
        @test dur(half_raster_rf, sys) ≈ delay(half_raster_rf, sys) +
            half_raster_shape_duration + sys.RF_ring_down_time

        # Explicit RF time shapes carry their nonzero first time in the stored Koma
        # delay; sys-aware timing recovers the Pulseq event delay and center.
        time_shaped_intervals = [3sys.RF_Δt, 4sys.RF_Δt]
        time_shape_start = 2sys.RF_Δt
        time_shape_center = time_shape_start + first(time_shaped_intervals)
        time_shaped_rf = RF(rf_samples, time_shaped_intervals, 0.0, rf_pulseq_delay + time_shape_start; center=time_shape_center - time_shape_start, use=Excitation())
        @test dwell(time_shaped_rf, sys) == time_shaped_intervals
        @test delay(time_shaped_rf, sys) ≈ rf_pulseq_delay
        @test rf_center(time_shaped_rf, sys) ≈ rf_center(time_shaped_rf) + delay(time_shaped_rf) - delay(time_shaped_rf, sys)
        time_shaped_duration = PulseDesigner.ceil_to_raster(
            dur(time_shaped_rf) - delay(time_shaped_rf, sys),
            sys.RF_Δt,
        )
        @test dur(time_shaped_rf, sys) ≈ delay(time_shaped_rf, sys) +
            time_shaped_duration + sys.RF_ring_down_time

        irregular_intervals = [3sys.RF_Δt / 2, 5sys.RF_Δt]
        irregular_time_shape_start = sys.RF_Δt / 2
        irregular_time_shape_center = 2sys.RF_Δt
        irregular_rf = RF(rf_samples, irregular_intervals, 0.0, rf_pulseq_delay + irregular_time_shape_start; center=irregular_time_shape_center - irregular_time_shape_start, use=Excitation())
        @test delay(irregular_rf, sys) ≈ rf_pulseq_delay
        @test rf_center(irregular_rf, sys) ≈ rf_center(irregular_rf) + delay(irregular_rf) - delay(irregular_rf, sys)
        irregular_shape_duration = PulseDesigner.ceil_to_raster(
            dur(irregular_rf) - delay(irregular_rf, sys),
            sys.RF_Δt,
        )
        @test dur(irregular_rf, sys) ≈ delay(irregular_rf, sys) +
            irregular_shape_duration + sys.RF_ring_down_time

        # A one-sample ADC still has a dwell interval in Pulseq.
        single_adc_dwell = sys.ADC_Δt
        single_adc_pulseq_delay = 10sys.ADC_Δt
        single_adc = ADC(1, single_adc_dwell, single_adc_pulseq_delay + single_adc_dwell / 2)
        @test dwell(single_adc, sys) ≈ single_adc_dwell
        @test delay(single_adc, sys) ≈ single_adc_pulseq_delay
        @test dur(single_adc, sys) ≈ single_adc_pulseq_delay + single_adc.N * single_adc_dwell + sys.ADC_dead_time
    end
    @testset "build_adc" begin
        sys = Scanner(ADC_Δt=2e-6, ADC_dead_time=6e-6, DUR_Δt=10e-6)
        # ADC builders accept Pulseq delay/duration and store Koma delay/sample span.
        adc_samples = 4
        adc_dwell = sys.ADC_Δt
        adc_pulseq_delay = 4e-6
        adc_duration = adc_samples * adc_dwell
        adc_freq_offset = 77.0
        adc_phase_offset = 0.2
        adc = PulseDesigner.make_adc(
            adc_samples, adc_dwell;
            delay=adc_pulseq_delay,
            freq_offset=adc_freq_offset,
            phase_offset=adc_phase_offset,
            sys,
        )
        @test adc.N == adc_samples
        @test adc.T == (adc_samples - 1) * adc_dwell
        @test adc.delay == max(adc_pulseq_delay, sys.ADC_dead_time) + adc_dwell / 2
        @test adc.Δf == adc_freq_offset
        @test adc.ϕ == adc_phase_offset

        seq = PulseDesigner.build_adc(adc_samples; duration=adc_duration, sys)
        adc = only(seq.ADC)
        @test adc.N == adc_samples
        @test adc.T == (adc_samples - 1) * adc_dwell
        @test adc.delay == sys.ADC_dead_time + adc_dwell / 2
        @test only(seq.DUR) == PulseDesigner.ceil_to_raster(dur(seq[1], sys), sys.DUR_Δt)
        single_dwell = 3e-6
        @test PulseDesigner.make_adc(1, single_dwell; sys).T == single_dwell

        @test_throws ErrorException PulseDesigner.make_adc(0, adc_dwell; sys)
        @test_throws ErrorException PulseDesigner.make_adc(adc_samples; sys)
        @test_throws ErrorException PulseDesigner.make_adc(
            adc_samples, adc_dwell; duration=adc_duration, sys,
        )
        @test_throws ErrorException PulseDesigner.make_adc(adc_samples, 0.0; sys)
    end
    @testset "build_arbitrary_rf" begin
        sys = Scanner(
            B1=Inf, Gmax=40u"mT/m", Smax=170u"T/m/s", RF_Δt=1u"μs",
            GR_Δt=10u"μs", RF_dead_time=0u"s", RF_ring_down_time=0u"s",
        )
        seq = PulseDesigner.build_arbitrary_rf(
            [1, 1, 1, 1], 90u"°"; dwell=10u"μs", bandwidth=2u"kHz",
            slice_thickness=5u"mm", sys,
        )
        gz, gzr = seq.GR.z
        slice_area = 2e3 * 40e-6 / (γ * 5e-3)
        @test gz.A * gz.T ≈ slice_area
        @test area(gzr) ≈ -slice_area * (1 - 20e-6 / 40e-6) -
            (area(gz) - slice_area) / 2

    end
    @testset "build_gauss_pulse" begin
        sys = Scanner(
            B1=Inf, Gmax=40u"mT/m", Smax=170u"T/m/s", RF_Δt=1u"μs",
            GR_Δt=10u"μs", RF_dead_time=0u"s", RF_ring_down_time=0u"s",
        )
        seq = PulseDesigner.build_gauss_pulse(
            90u"°"; duration=1u"ms", bandwidth=2u"kHz",
            slice_thickness=5u"mm", max_grad=30u"mT/m", max_slew=120u"T/m/s", sys,
        )
        gz, gzr = seq.GR.z
        slice_area = 2e3 * 1e-3 / (γ * 5e-3)
        @test gz.A * gz.T ≈ slice_area
        @test area(gzr) ≈ -slice_area / 2 - (area(gz) - slice_area) / 2
    end
    @testset "build_adiabatic_pulse" begin
        sys = Scanner(
            B1=Inf, Gmax=40u"mT/m", Smax=170u"T/m/s", RF_Δt=1u"μs",
            GR_Δt=10u"μs", RF_dead_time=0u"s", RF_ring_down_time=0u"s",
        )
        seq = PulseDesigner.build_adiabatic_pulse(
            :wurst; duration=4u"ms", dwell=2u"μs", bandwidth=6u"kHz",
            freq_offset=100u"Hz", phase_offset=30u"°", sys,
        )
        rf = seq.RF[1, 1]
        @test rf.Δf == 100.0
        @test rf.ϕ ≈ π / 6 atol=1e-5

        seq = PulseDesigner.build_adiabatic_pulse(
            :wurst; duration=4u"ms", bandwidth=6u"kHz", slice_thickness=5u"mm", sys,
        )
        rf = seq.RF[1, 1]
        gz, gzr = seq.GR.z
        slice_area = 6e3 * 4e-3 / (γ * 5e-3)
        center_pos = rf_center(rf, sys) / 4e-3
        @test gz.A * gz.T ≈ slice_area
        @test area(gzr) ≈ -slice_area * (1 - center_pos) - (area(gz) - slice_area) / 2
    end
    @testset "build_sinc_pulse" begin
        sys = Scanner(
            B1=Inf, Gmax=40u"mT/m", Smax=170u"T/m/s", RF_Δt=1u"μs",
            GR_Δt=10u"μs", RF_dead_time=0u"s", RF_ring_down_time=0u"s",
        )
        flip_angle = 90u"°"
        duration = 1u"ms"
        seq = PulseDesigner.build_sinc_pulse(
            flip_angle; duration, slice_thickness=5u"mm", freq_offset=1u"kHz", sys,
        )
        rf = seq.RF[1, 1]
        gz, gzr = seq.GR.z
        slice_area = 4 / (γ * 5e-3)
        @test gz.A * gz.T ≈ slice_area
        @test area(gzr) ≈ -slice_area / 2 - (area(gz) - slice_area) / 2
        @test rf.delay ≈ gz.rise + gz.delay + sys.RF_Δt / 2
        @test rf.Δf == 1e3

    end
    @testset "Spiral" begin
        sys = Scanner()
        sys.Smax = 150    # [mT/m/ms]
        sys.Gmax = 500e-3 # [T/m]
        sys.GR_Δt = 4e-6  # [s]
        FOV = 0.2       # [m]
        N = 80          # Reconstructed image N×N
        Nint = 8
        λ = 2.1
        spiral = PulseDesigner.spiral_base(FOV, N, sys; λ=λ, BW=120e3, Nint)
        # Look at the k_space generated
        @test spiral(0).DEF["λ"] ≈ λ
    end
    @testset "Radial" begin
        sys = Scanner()
        N = 80
        Nspokes = ceil(Int64, π/2 * N ) #Nyquist in the radial direction
        FOV = 0.2
        spoke = PulseDesigner.radial_base(FOV, N, sys)
        @test spoke.DEF["Δθ"] ≈ π / Nspokes
    end
end

@testitem "Motion" tags=[:base] begin
    @testset "Constructors" begin
        action = Rotate(10.0, 20.0, 40.0, (0.0, 0.0, 0.0))
        spins = AllSpins()
        # TimeCurve constructors
        time = TimeRange(t_start=0.0, t_end=1.0)
        time = Periodic(period=1.0, asymmetry=0.5)
        time = Periodic(period=1.0, asymmetry=0.0)
        time = Periodic(period=1.0, asymmetry=1.0)
        time = TimeCurve([0.0, eps()], [0.0, 1.0])
        # Motion constructors
        m = Motion(action, time, spins)
        @test Motion(action) == m
        @test Motion(action, time) == m
        @test Motion(action, spins) == m
        # MotionList
        @test MotionList() == NoMotion()
        @test MotionList(m) == m
    end
    @testset "Times" begin
        tr = translate(0.0, 0.1, 0.2, TimeRange(0.0, 1.0))
        rt = rotate(10.0, 20.0, 30.0, TimeCurve(t=[0.0, 0.5, 0.8], t_unit=[0.0, 0.5, 1.0]))
        ml = MotionList(rt, tr)
        @test times(tr) == [0.0, 1.0]
        @test times(rt) == [0.0, 0.5, 0.8]
        @test times(ml) == [0.0, 0.5, 0.8, 1.0]
    end
    @testset "Subset" begin
        rt = Rotate(10.0, 20.0, 40.0, (0.0, 0.0, 0.0))
        tr = Translate(0.1, 0.2, 0.3)
        time = TimeRange(0.0, eps())
        spins = AllSpins()
        rng = rng = 1:2:5
        # NoMotion
        nm = NoMotion()
        @test nm[rng] == nm
        # Motion
        m = Motion(rt, time, spins)
        @test m[rng] == m
        # MotionList 
        ml = MotionList(Motion(rt, time, spins), Motion(tr, time, spins))
        @test ml[rng] == ml 
    end
    # Spin Positions
    @testset "NoMotion" begin
        ph = Phantom(x=[1.0, 2.0], y=[1.0, 2.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        xt, yt, zt = get_spin_coords(ph.motion, ph.x, ph.y, ph.z, t')
        @test xt == ph.x
        @test yt == ph.y
        @test zt == ph.z
    end
    @testset "Translate" begin
        ph = Phantom(x=[1.0], y=[1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        dx, dy, dz = [1.0, 0.0, 0.0]
        vx, vy, vz = [dx, dy, dz] ./ (t_end - t_start)
        translation = translate(dx, dy, dz, TimeRange(t_start, t_end))
        xt, yt, zt = get_spin_coords(translation, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ vx.*t'
        @test yt == ph.y .+ vy.*t'
        @test zt == ph.z .+ vz.*t'
    end
    @testset "PeriodicTranslate" begin
        ph = Phantom(x=[1.0], y=[1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        period = 2.0
        asymmetry = 0.5
        dx, dy, dz = [1.0, 0.0, 0.0]
        vx, vy, vz = [dx, dy, dz] ./ (t_end - t_start)
        periodictranslation = translate(dx, dy, dz, Periodic(period=period, asymmetry=asymmetry))
        xt, yt, zt = get_spin_coords(periodictranslation, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ vx.*t'
        @test yt == ph.y .+ vy.*t'
        @test zt == ph.z .+ vz.*t'
    end
    @testset "Rotate" begin
        # Simple-axis Rotation Constructors
        @test RotateX(90.0) == Rotate(90.0, 0.0, 0.0, CenterOfMass())
        @test RotateY(90.0) == Rotate(0.0, 90.0, 0.0, CenterOfMass())
        @test RotateZ(90.0) == Rotate(0.0, 0.0, 90.0, CenterOfMass())
        # Test get_spin_coords
        ph = Phantom(x=[1.0, 1.0, -1.0, -1.0], y=[1.0, -1.0, 1.0, -1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        pitch, roll, yaw = 45.0, 45.0, 45.0
        # One single rotation (around center of mass)
        rotation = rotate(pitch, roll, yaw, TimeRange(t_start, t_end))
        xt, yt, zt = get_spin_coords(rotation, ph.x, ph.y, ph.z, t')
        R = rotz(π*yaw/180) * roty(π*roll/180) * rotx(π*pitch/180)
        r = hcat(ph.x, ph.y, ph.z)'
        rotated = R * r 
        rot_x, rot_y, rot_z = eachrow(rotated)
        @test xt[: ,end] ≈ rot_x
        @test yt[: ,end] ≈ rot_y
        @test zt[: ,end] ≈ rot_z
        # One single rotation (around displaced center)
        center = (0.1, 0.2, 0.3)
        rotation_displaced = rotate(pitch, roll, yaw, TimeRange(t_start, t_end); center=center)
        @test !(rotation ≈ rotation_displaced) & !(rotation_displaced ≈ rotation)
        xt, yt, zt = get_spin_coords(rotation_displaced, ph.x, ph.y, ph.z, t')
        R = rotz(π*yaw/180) * roty(π*roll/180) * rotx(π*pitch/180)
        r = hcat(ph.x .- center[1], ph.y .- center[2], ph.z .- center[3])'
        rotated = R * r 
        rot_x, rot_y, rot_z = eachrow(rotated)
        @test xt[: ,end] ≈ rot_x .+ center[1]
        @test yt[: ,end] ≈ rot_y .+ center[2]
        @test zt[: ,end] ≈ rot_z .+ center[3]
        # Check if two consecutive rotations (α and β) produce the same result as a single (α + β) rotation
        t = [1.0] 
        r1 = MotionList(
            rotate(0.0, 0.0, yaw/2, TimeRange(t_start, t_end/2)),
            rotate(0.0, 0.0, yaw/2, TimeRange(t_end/2, t_end))
        )
        r2 = rotate(0.0, 0.0, yaw, TimeRange(t_start, t_end))
        xt1, yt1, zt1 = get_spin_coords(r1, ph.x, ph.y, ph.z, t)
        xt2, yt2, zt2 = get_spin_coords(r2, ph.x, ph.y, ph.z, t)
        @test xt1 ≈ xt2
        @test yt1 ≈ yt2
        @test zt1 ≈ zt2
    end
    @testset "PeriodicRotation" begin
        ph = Phantom(x=[1.0, 1.0, -1.0, -1.0], y=[1.0, -1.0, 1.0, -1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        period = 2.0
        asymmetry = 0.5
        pitch = 45.0
        roll = 45.0
        yaw = 45.0
        periodicrotation = rotate(pitch, roll, yaw, Periodic(period=period, asymmetry=asymmetry))
        xt, yt, zt = get_spin_coords(periodicrotation, ph.x, ph.y, ph.z, t')
        R = rotz(π*yaw/180) * roty(π*roll/180) * rotx(π*pitch/180)
        r = hcat(ph.x, ph.y, ph.z)'
        rotated = R * r 
        rot_x, rot_y, rot_z = eachrow(rotated)
        @test xt[: ,end] ≈ rot_x
        @test yt[: ,end] ≈ rot_y
        @test zt[: ,end] ≈ rot_z
    end
    @testset "HeartBeat" begin
        ph = Phantom(x=[1.0], y=[1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        circumferential_strain = -0.1
        radial_strain = 0.0
        longitudinal_strain = -0.1
        hb = heartbeat(circumferential_strain, radial_strain, longitudinal_strain, TimeRange(t_start, t_end))
        xt, yt, zt = get_spin_coords(hb, ph.x, ph.y, ph.z, t')
        r = sqrt.(ph.x .^ 2 + ph.y .^ 2)
        θ = atan.(ph.y, ph.x)
        @test xt[:,end] == ph.x .* (1 .+ circumferential_strain * maximum(r) .* cos.(θ))
        @test yt[:,end] == ph.y .* (1 .+ circumferential_strain * maximum(r) .* sin.(θ))
        @test zt[:,end] == ph.z .* (1 .+ longitudinal_strain)
    end
    @testset "PeriodicHeartBeat" begin
        ph = Phantom(x=[1.0], y=[1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        period = 2.0
        asymmetry = 0.5
        circumferential_strain = -0.1
        radial_strain = 0.0
        longitudinal_strain = -0.1
        periodic_hb = heartbeat(circumferential_strain, radial_strain, longitudinal_strain, Periodic(period=period, asymmetry=asymmetry))
        xt, yt, zt = get_spin_coords(periodic_hb, ph.x, ph.y, ph.z, t')
        r = sqrt.(ph.x .^ 2 + ph.y .^ 2)
        θ = atan.(ph.y, ph.x)
        @test xt[:,end] == ph.x .* (1 .+ circumferential_strain * maximum(r) .* cos.(θ))
        @test yt[:,end] == ph.y .* (1 .+ circumferential_strain * maximum(r) .* sin.(θ))
        @test zt[:,end] == ph.z .* (1 .+ longitudinal_strain)
    end
    @testset "Path" begin
        # 1 spin
        ph = Phantom(x=[1.0], y=[1.0])
        Ns = length(ph)
        t_start = 0.0
        t_end = 1.0
        Nt = 10
        dx = rand(Ns, Nt)
        dy = rand(Ns, Nt)
        dz = rand(Ns, Nt)
        pt = path(dx, dy, dz, TimeRange(t_start, t_end))
        t = range(t_start, t_end, Nt)
        xt, yt, zt = get_spin_coords(pt, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ dx
        @test yt == ph.y .+ dy
        @test zt == ph.z .+ dz
        # More than 1 spin
        ph = Phantom(x=[1.0, 2.0], y=[1.0, 2.0])
        Ns = length(ph)
        t_start = 0.0
        t_end = 1.0
        Nt = 10
        dx = rand(Ns, Nt)
        dy = rand(Ns, Nt)
        dz = rand(Ns, Nt)
        pt = path(dx, dy, dz, TimeRange(t_start, t_end))
        t = range(t_start, t_end, Nt)
        xt, yt, zt = get_spin_coords(pt, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ dx
        @test yt == ph.y .+ dy
        @test zt == ph.z .+ dz
    end
    @testset "FlowPath" begin
        # 1 spin
        ph = Phantom(x=[1.0], y=[1.0])
        Ns = length(ph)
        t_start = 0.0
        t_end = 1.0
        Nt = 10
        dx = rand(Ns, Nt)
        dy = rand(Ns, Nt)
        dz = rand(Ns, Nt)
        fp = flowpath(dx, dy, dz, Bool.(zeros(Ns, Nt)), TimeRange(t_start, t_end))
        t = range(t_start, t_end, Nt)
        xt, yt, zt = get_spin_coords(fp, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ dx
        @test yt == ph.y .+ dy
        @test zt == ph.z .+ dz
        # More than 1 spin
        ph = Phantom(x=[1.0, 2.0], y=[1.0, 2.0])
        Ns = length(ph)
        t_start = 0.0
        t_end = 1.0
        Nt = 10
        dx = rand(Ns, Nt)
        dy = rand(Ns, Nt)
        dz = rand(Ns, Nt)
        fp = flowpath(dx, dy, dz, Bool.(zeros(Ns, Nt)), TimeRange(t_start, t_end))
        t = range(t_start, t_end, Nt)
        xt, yt, zt = get_spin_coords(fp, ph.x, ph.y, ph.z, t')
        @test xt == ph.x .+ dx
        @test yt == ph.y .+ dy
        @test zt == ph.z .+ dz
    end
    @testset "Translate + Rotate" begin
        ph = Phantom(x=[1.0, 1.0, -1.0, -1.0], y=[1.0, -1.0, 1.0, -1.0])
        t_start=0.0; t_end=1.0 
        t = collect(range(t_start, t_end, 11))
        # Translate
        dx, dy, dz = [1.0, 0.0, 0.0]
        vx, vy, vz = [dx, dy, dz] ./ (t_end - t_start)
        translation = translate(dx, dy, dz, TimeRange(t_start, t_end))
        # Rotate
        pitch, roll, yaw = [45.0, 45.0, 45.0]
        rotation = rotate(pitch, roll, yaw, TimeRange(t_start, t_end))
        R = rotz(π*yaw/180) * roty(π*roll/180) * rotx(π*pitch/180)
        r = hcat(ph.x, ph.y, ph.z)'
        rotated = R * r 
        rot_x, rot_y, rot_z = eachrow(rotated)
        # Combination into a MotionList
        motion = MotionList(translation, rotation)
        xt, yt, zt = get_spin_coords(motion, ph.x, ph.y, ph.z, t')
        @test xt[: ,end] ≈ rot_x .+ vx*t[end]
        @test yt[: ,end] ≈ rot_y .+ vy*t[end]
        @test zt[: ,end] ≈ rot_z .+ vz*t[end]
    end
    @testset "Key Time Points" begin
        # Sequence duration
        t_start = 0.0
        t_end   = 1.5 
        # TimeCurve parameters
        t        = [0.0, 0.1, 0.3]
        t_unit   = [0.0, 0.4, 1.0]
        periods  = [1.0, 0.5, 2.0]
        dx = dy  = [0.0 0.0 0.0 0.0]
        dz       = [3.0 4.0 -4. -3.]
        reset    = [false false true false]
        ϵ = KomaMRIBase.MIN_RISE_TIME

        # Key time points ("manually" determined):
        # Periodic case
        period_times_p  = [t+δ for t in (0.0, 0.3, 0.45, 1.05, 1.35, 1.5) for δ in (-ϵ, ϵ) if (t+δ) > t_start && (t+δ) < t_end]
        reset_times_p   = [0.2, 0.4, 0.85, 1.25, 1.45] .- ϵ
        # Non-periodic case:
        period_times_np = [t+δ for t in (0.0, 0.3, 0.45, 1.05) for δ in (-ϵ, ϵ) if (t+δ) > t_start && (t+δ) < 1.05]
        reset_times_np  = [0.2, 0.4, 0.85] .- ϵ

        # Any motion with no spin resets (only key time points derived from periods):
        pth  = path(dx, dy, dz, TimeCurve(t, t_unit, true, periods), AllSpins())
        seqd_t = [t_start, t_end]
        KomaMRIBase.add_key_time_points!(seqd_t, pth)
        @test seqd_t ≈ [t_start; t_end; period_times_p]
        # FlowPath with a spin reset (key time points derived from both periods and spin resets):
        fpth = flowpath(dx, dy, dz, reset, TimeCurve(t, t_unit, true, periods), AllSpins())
        seqd_t = [t_start, t_end]
        KomaMRIBase.add_key_time_points!(seqd_t, fpth)
        @test sort(seqd_t) ≈ sort([t_start; t_end; period_times_p; reset_times_p])

        # MotionList 
        # (periodic case)
        ml = MotionList(pth, fpth)
        seqd_t = [t_start, t_end]
        KomaMRIBase.add_key_time_points!(seqd_t, ml)
        @test sort(unique(seqd_t)) ≈ sort([t_start; t_end; period_times_p; reset_times_p])
        # (non-periodic case)
        pth  = path(dx, dy, dz, TimeCurve(t, t_unit, false, periods), AllSpins())
        fpth = flowpath(dx, dy, dz, reset, TimeCurve(t, t_unit, false, periods), AllSpins())
        ml = MotionList(pth, fpth)
        seqd_t = [t_start, t_end]
        KomaMRIBase.add_key_time_points!(seqd_t, ml)
        @test unique(seqd_t) ≈ [t_start; t_end; period_times_np; reset_times_np]
    end
end

@testitem "Phantom" tags = [:base] begin
    using Suppressor
    # Phantom Struct Fields
    name = "Bulks"
    x = [-2e-3; -1e-3; 0.0; 1e-3; 2e-3]
    y = [-4e-3; -2e-3; 0.0; 2e-3; 4e-3]
    z = [-6e-3; -3e-3; 0.0; 3e-3; 6e-3]
    ρ = [0.2; 0.4; 0.6; 0.8; 1.0]
    T1 = [0.9; 0.9; 0.5; 0.25; 0.4]
    T2 = [0.09; 0.05; 0.04; 0.07; 0.005]
    T2s = [0.1; 0.06; 0.05; 0.08; 0.015]
    Δw = [-2e-6; -1e-6; 0.0; 1e-6; 2e-6]
    Dλ1 = [-4e-6; -2e-6; 0.0; 2e-6; 4e-6]
    Dλ2 = [-6e-6; -3e-6; 0.0; 3e-6; 6e-6]
    Dθ = [-8e-6; -4e-6; 0.0; 4e-6; 8e-6]
    # Motion
    Ns = length(x)
    Nt = 3
    t_start = 0.0
    t_end = 1.0
    tr = translate(0.05, 0.05, 0.0, Periodic(period=0.5, asymmetry=0.5))
    rt = rotate(0.0, 0.0, 90.0, TimeRange(t_start=0.05, t_end=0.5), SpinRange(1:3))
    pt = path(0.01 .* rand(Ns, Nt), 0.01 .* rand(Ns, Nt), 0.01 .* rand(Ns, Nt), TimeRange(t_start, t_end), SpinRange(2:2:4))
    @testset "Comparison" begin
        obj1 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, motion=MotionList(tr, rt))
        obj2 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, motion=MotionList(tr, rt))
        @test obj1 == obj2
        obj2.x .+= 1e-10
        @test obj1 ≈ obj2
        obj1.motion = NoMotion()
        @test !(obj1 == obj2)
        @test !(obj1  ≈ obj2)
    end
    @testset "Size and Length" begin
        obj1 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, motion=MotionList(tr, rt))
        @test size(obj1) == size(ρ)
        @test length(obj1) == length(ρ)
    end
    @testset "Subset" begin 
        motion = MotionList(tr, rt)
        obj1 = Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, motion)
        rng = 1:2:5
        obj2 = Phantom(
            name, x[rng], y[rng], z[rng], 
            ρ[rng], T1[rng], T2[rng], T2s[rng], 
            Δw[rng], Dλ1[rng], Dλ2[rng], Dθ[rng], 
            motion[rng]
        )
        # Phantom subset
        @test obj1[rng] == obj2
        @test @view(obj1[rng]) == obj2
        # Phantom view
        obj_view = @view(obj1[rng])
        obj_view.ρ .= 0.0
        @test obj_view.ρ == obj1[rng].ρ
        # BitVector range
        obj3 = copy(obj1)
        rng = obj1.x .> 0
        obj1.motion = translate(5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), SpinRange(rng))
        obj3.motion = translate(5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), SpinRange(1:length(obj3)))
        @test obj1[rng] == obj3[rng]
        @test obj1[rng].motion == obj3.motion[rng]
    end
    @testset "Addition" begin
        obj1 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ)
        rng = 1:2:5
        obj2 = obj1[rng]
        oba = Phantom(
            name, [x; x[rng]], [y; y[rng]], [z; z[rng]], 
            [ρ; ρ[rng]], [T1; T1[rng]], [T2; T2[rng]], [T2s; T2s[rng]], 
            [Δw; Δw[rng]], [Dλ1; Dλ1[rng]], [Dλ2; Dλ2[rng]], [Dθ; Dθ[rng]], 
            vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        )
        # NOTE: these vcat methods must be simplified once the Vector{<:Motion} approach is accomplished: 
        # https://github.com/JuliaHealth/KomaMRI.jl/issues/480
        # NoMotion + NoMotion
        @test obj1 + obj2 == oba
        # NoMotion + MotionList
        obj2.motion = MotionList(tr, rt)
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # MotionList + NoMotion
        obj1.motion = MotionList(tr, rt)
        obj2.motion = NoMotion()
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # NoMotion + Motion
        obj1.motion = NoMotion()
        obj2.motion = tr
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # Motion + NoMotion
        obj1.motion = tr
        obj2.motion = NoMotion()
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # MotionList + MotionList
        obj1.motion = MotionList(tr, rt)
        obj2.motion = MotionList(tr, rt)
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # Motion + Motion
        obj1.motion = tr
        obj2.motion = rt
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # Motion + MotionList
        obj1.motion = tr
        obj2.motion = MotionList(tr, rt)
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
        # MotionList + Motion
        obj1.motion = MotionList(tr, rt)
        obj2.motion = tr
        oba.motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2))
        @test obj1 + obj2 == oba
    end 
    @testset "Scalar multiplication" begin
        obj1 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, motion=MotionList(tr, rt))
        c = 7
        obc = Phantom(name=name, x=x, y=y, z=z, ρ=c*ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, motion=MotionList(tr, rt))
        @test c * obj1 == obc
    end
    @testset "Brain Phantom 2D" begin
        ph = brain_phantom2D()
        @test ph.name == "brain2D_axial"
        @test KomaMRIBase.get_dims(ph) == Bool[1, 1, 0]
    end
    @testset "Brain Phantom 3D" begin
        ph = brain_phantom3D()
        @test ph.name == "brain3D"
        @test KomaMRIBase.get_dims(ph) == Bool[1, 1, 1]
    end
    @testset "Pelvis Phantom" begin
        ph = pelvis_phantom2D()
        @test ph.name == "pelvis2D"
        @test KomaMRIBase.get_dims(ph) == Bool[1, 1, 0]
    end
    @testset "Heart Phantom" begin
        ph = heart_phantom()
        @test ph.name == "LeftVentricle"
        @test KomaMRIBase.get_dims(ph) == Bool[1, 1, 0]
        @test 0 < length(heart_phantom(spins_per_voxel=20)) < length(ph)
    end
end

@testitem "Scanner" tags=[:base] begin
    B0, B1, Gmax, Smax = 1.5, 10e-6, 60e-3, 500
    ADC_Δt, DUR_Δt, GR_Δt, RF_Δt = 2e-6, 1e-5, 1e-5, 1e-6
    RF_ring_down_time, RF_dead_time, ADC_dead_time = 20e-6, 100e-6, 10e-6
    sys = Scanner(B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt, RF_ring_down_time, RF_dead_time, ADC_dead_time)
    @test sys.B0 ≈ B0 && sys.B1 ≈ B1 && sys.Gmax ≈ Gmax && sys.Smax ≈ Smax
end

@testitem "TrapezoidalIntegration" tags=[:base] begin
    dt = Float64[1 1 1 1]
    x  = Float64[0 1 2 1 0]
    @test KomaMRIBase.trapz(dt, x)[1] ≈ 4 #Triangle area = bh/2, with b = 4 and h = 2
    @test KomaMRIBase.cumtrapz(dt, x) ≈ [0.5 2 3.5 4]
end
