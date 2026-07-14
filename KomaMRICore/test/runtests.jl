using TestItems, TestItemRunner

### NOTE: by default, tests are run on the CPU with the number of threads set to
#   Threads.nthreads(). To run on a specific GPU backend, add the name of the
#   backend package ("AMDGPU", "CUDA", "Metal", or "oneAPI") to the test/Project.toml
#   file in KomaMRICore and pass the name as a test argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; test_args=["CUDA"])
#
#   To run on the cpu with a specific number of threads, pass the number of threads
#   as a julia argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; julia_args=`--threads=4`)
#
#   For changing the default backend used for testing,
#   modify the [preferences.KomaMRICore] section in the test/Project.toml:
#
#     [preferences.KomaMRICore]
#     test_backend = "CPU"
#
#   For the backend preference to take effect, you need to:
#   - REPL testing: No action needed. `] test` should pick up the preference right away.
#   - VSCode testing: You need to restart VSCode.
#
#   Sadly, LocalPreferences.toml are not picked up by VScode (that could be .gitignore'd),
#   so we had put them into the test/Project.toml.
#
###

#Environment variable set by CI
const CI    = get(ENV, "CI", nothing)
const group = get(ENV, "TEST_GROUP", :core) |> Symbol

@run_package_tests filter=ti->(group in ti.tags)&&(isnothing(CI) || :skipci ∉ ti.tags) #verbose=true

@testitem "simulation-method sampling rule default" tags=[:core, :nomotion] begin
    struct RuleDefaultMethod <: KomaMRICore.SimulationMethod end
    struct RuleDefaultTestRule <: SamplingRule end

    KomaMRICore.default_sampling_rule(::RuleDefaultMethod, sim_params) = RuleDefaultTestRule()

    sim_params = KomaMRICore.default_sim_params(Dict{String,Any}("sim_method" => RuleDefaultMethod()))
    @test !haskey(sim_params, "sampling_rule")
    rule = KomaMRICore.simulation_sampling_rule(sim_params["sim_method"], sim_params)
    @test rule isa KomaMRICore.IntegrationNodeSamplingRule
    @test rule.rule isa RuleDefaultTestRule
    @test KomaMRICore.integration_nodes(BlochMagnusConst1()) == (0,)
    @test KomaMRICore.integration_nodes(BlochMagnusLin2()) == (0, 1)
    @test KomaMRICore.integration_nodes(BlochMagnusMid2()) == (1//2,)
    @test KomaMRICore.integration_nodes(BlochMagnusLinComm2()) == (0, 1)
    @test KomaMRICore.integration_nodes(BlochMagnusQuad4()) == (0, 1//2, 1)
    @test KomaMRICore.eval_intervals_per_step(BlochMagnusConst1()) == 1
    @test KomaMRICore.eval_intervals_per_step(BlochMagnusMid2()) == 2
    @test KomaMRICore.eval_intervals_per_step(BlochMagnusQuad4()) == 2
    @test KomaMRICore.eval_intervals_per_step(BlochMagnusGL4()) == 3

    override = MaxStepSizeRule(2e-3, 3e-5)
    sim_params = KomaMRICore.default_sim_params(Dict{String,Any}(
        "sim_method" => RuleDefaultMethod(),
        "sampling_rule" => override,
    ))
    @test sim_params["sampling_rule"] === override
    rule = KomaMRICore.simulation_sampling_rule(sim_params["sim_method"], sim_params)
    @test rule.rule === override

    sim_params = KomaMRICore.default_sim_params(Dict{String,Any}("sim_method" => BlochMagnusQuad4()))
    sim_params["Δt_rf"] = 2.5e-4
    rule = KomaMRICore.simulation_sampling_rule(sim_params["sim_method"], sim_params)
    @test rule.rule.Δt_rf == 2.5e-4

    sim_params["preserve_samples"] = ()
    rule = KomaMRICore.simulation_sampling_rule(sim_params["sim_method"], sim_params)
    @test KomaMRIBase.preserved_samples(rule) == ()

    delays = Sequence()
    delays += Delay(1e-3)
    delays += Delay(1e-3)
    base_seqd = discretize(delays; sampling_rule=MaxStepSizeRule(Inf, Inf))
    @test base_seqd.Δt ≈ [2e-3]

    method_rule = KomaMRICore.simulation_sampling_rule(BlochMagnusQuad4(), Dict("Δt" => Inf, "Δt_rf" => Inf))
    method_seqd = discretize(delays; sampling_rule=method_rule)
    @test method_seqd.t ≈ [0.0, 2e-3]
    @test method_seqd.Δt ≈ [2e-3]
end

@testitem "BlochMagnusQuad4 midpoint sampling" tags=[:core, :nomotion] begin
    # Quadratic RF intervals are represented as start, true midpoint, and stop rows.
    seq = Sequence()
    @addblock seq += RF([1.0, 2.0, 4.0] .* 1e-6, 1e-3, [0.0, 500.0, 1000.0])
    rule = KomaMRICore.simulation_sampling_rule(BlochMagnusQuad4(), Dict("Δt" => Inf, "Δt_rf" => 2.5e-4))
    seqd = discretize(seq; motion=NoMotion(), sampling_rule=rule)

    let i = firstindex(seqd.Δt)
        while i <= lastindex(seqd.Δt)
            if seqd.excitation_bool[i] && seqd.Δt[i] > 0
                midpoint = i + 1
                stop = i + 2
                @test midpoint <= lastindex(seqd.Δt)
                @test stop <= lastindex(seqd.t)
                @test seqd.excitation_bool[midpoint]
                @test seqd.Δt[i] ≈ seqd.Δt[midpoint]
                @test seqd.t[midpoint] ≈ (seqd.t[i] + seqd.t[stop]) / 2

                direct = KomaMRIBase.evaluate_sequence_at(seq, [seqd.t[midpoint]])
                @test only(direct.B1) ≈ seqd.B1[midpoint]
                @test only(direct.Δf) ≈ seqd.Δf[midpoint]
                i += 2
            else
                i += 1
            end
        end
    end
end

@testitem "BlochMagnusGL RF sampling" tags=[:core, :nomotion] begin
    # GL RF intervals keep boundary rows but use two interior Gauss-Legendre nodes.
    seq = Sequence()
    @addblock seq += RF([1.0, 2.0, 4.0] .* 1e-6, 1e-3, [0.0, 500.0, 1000.0])
    rule = KomaMRICore.simulation_sampling_rule(BlochMagnusGL4(), Dict("Δt" => Inf, "Δt_rf" => 2.5e-4))
    seqd = discretize(seq; motion=NoMotion(), sampling_rule=rule)
    c_minus = 1 / 2 - sqrt(3) / 6
    c_plus = 1 / 2 + sqrt(3) / 6

    let i = firstindex(seqd.Δt)
        while i <= lastindex(seqd.Δt)
            if seqd.excitation_bool[i] && seqd.Δt[i] > 0
                i_minus, i_plus, stop = i + 1, i + 2, i + 3
                @test i_plus <= lastindex(seqd.Δt)
                @test stop <= lastindex(seqd.t)
                @test seqd.excitation_bool[i_minus]
                @test seqd.excitation_bool[i_plus]
                Δt = sum(seqd.Δt[i:i_plus])
                @test seqd.t[i_minus] ≈ seqd.t[i] + c_minus * Δt
                @test seqd.t[i_plus] ≈ seqd.t[i] + c_plus * Δt

                direct = KomaMRIBase.evaluate_sequence_at(seq, [seqd.t[i_minus], seqd.t[i_plus]])
                @test direct.B1 ≈ seqd.B1[[i_minus, i_plus]]
                @test direct.Δf ≈ seqd.Δf[[i_minus, i_plus]]
                i += 3
            else
                i += 1
            end
        end
    end
end

@testitem "BlochMagnus RF block length counts simulation steps" tags=[:core, :nomotion] begin
    seq = Sequence()
    @addblock seq += RF([1.0, 2.0, 4.0] .* 1e-6, 1e-3, [0.0, 500.0, 1000.0])

    for (sim_method, eval_stride) in ((BlochMagnusMid2(), 2), (BlochMagnusQuad4(), 2), (BlochMagnusGL4(), 3))
        rule = KomaMRICore.simulation_sampling_rule(sim_method, Dict("Δt" => Inf, "Δt_rf" => 2.5e-4))
        seqd = discretize(seq; motion=NoMotion(), sampling_rule=rule)
        ranges, is_excitation = KomaMRICore.get_sim_ranges(
            seqd;
            max_rf_block_length=1,
            eval_intervals_per_step=KomaMRICore.eval_intervals_per_step(sim_method),
        )
        @test all(length(r) - 1 == eval_stride for (r, excitation) in zip(ranges, is_excitation) if excitation)
    end
end

@testitem "BlochMagnus simulates RF-center/max-step collision" tags=[:core, :nomotion] begin
    sys = Scanner()
    sys.limits.Smax = 100.0
    rf_B1 = 4.9e-6
    rf_duration = 3.2e-3
    slice_thickness = 8e-3
    slice_bandwidth = 5e3
    Gz = slice_bandwidth / (γ * slice_thickness)
    seq = PulseDesigner.RF_sinc(rf_B1, rf_duration, sys; G=[Gz, 0.0, 0.0], TBP=8)
    obj = Phantom(x=[0.0], y=[0.0], z=[0.0], ρ=[1.0], T1=[Inf], T2=[Inf])

    for sim_method in (BlochMagnusMid2(), BlochMagnusQuad4(), BlochMagnusGL2(), BlochMagnusGL4())
        sim_params = KomaMRICore.default_sim_params(Dict{String,Any}(
            "gpu" => false,
            "precision" => "f64",
            "return_type" => "state",
            "sim_method" => sim_method,
            "Δt" => 1e-3,
            "Δt_rf" => 100e-6,
        ))
        simulate(obj, seq, sys; sim_params, verbose=false)
        @test true
    end
end

@testitem "Spinors×Mag" tags=[:core, :nomotion] begin
    using KomaMRICore: Rx, Ry, Rz, Q, rotx, roty, rotz, Un, Rφ, Rg

    ## Verifying that operators perform counter-clockwise rotations
    v = [1, 2, 3]
    m = Mag([complex(v[1:2]...)], [v[3]])
    # Rx
    @test rotx(π/2) * v    ≈ [1, -3, 2]
    @test (Rx(π/2) * m).xy ≈ [1.0 - 3.0im]
    @test (Rx(π/2) * m).z  ≈ [2.0]
    # Ry
    @test roty(π/2) * v    ≈ [3, 2, -1]
    @test (Ry(π/2) * m).xy ≈ [3.0 + 2.0im]
    @test (Ry(π/2) * m).z  ≈ [-1.0]
    # Rz
    @test rotz(π/2) * v    ≈ [-2, 1, 3]
    @test (Rz(π/2) * m).xy ≈ [-2.0 + 1.0im]
    @test (Rz(π/2) * m).z  ≈ [3.0]
    # Rn
    @test Un(π/2, [1,0,0]) * v ≈ rotx(π/2) * v
    @test Un(π/2, [0,1,0]) * v ≈ roty(π/2) * v
    @test Un(π/2, [0,0,1]) * v ≈ rotz(π/2) * v
    @test (Q(π/2, 1.0+0.0im, 0.0) * m).xy ≈ (Rx(π/2) * m).xy
    @test (Q(π/2, 1.0+0.0im, 0.0) * m).z  ≈ (Rx(π/2) * m).z
    @test (Q(π/2, 0.0+1.0im, 0.0) * m).xy ≈ (Ry(π/2) * m).xy
    @test (Q(π/2, 0.0+1.0im, 0.0) * m).z  ≈ (Ry(π/2) * m).z
    @test (Q(π/2, 0.0+0.0im, 1.0) * m).xy ≈ (Rz(π/2) * m).xy
    @test (Q(π/2, 0.0+0.0im, 1.0) * m).z  ≈ (Rz(π/2) * m).z

    ## Verify that Spinor rotation = matrix rotation
    v = rand(3)
    n = rand(3); n = n ./ sqrt(sum(n.^2))
    m = Mag([complex(v[1:2]...)], [v[3]])
    φ, θ, φ1, φ2 = rand(4) * 2π
    # Rx
    vx = rotx(θ) * v
    mx = Rx(θ) * m
    @test [real(mx.xy); imag(mx.xy); mx.z] ≈ vx
    # Ry
    vy = roty(θ) * v
    my = Ry(θ) * m
    @test [real(my.xy); imag(my.xy); my.z] ≈ vy
    # Rz
    vz = rotz(θ) * v
    mz = Rz(θ) * m
    @test [real(mz.xy); imag(mz.xy); mz.z] ≈ vz
    # Rφ
    vφ = Un(θ, [sin(φ); cos(φ); 0.0]) * v
    mφ = Rφ(φ,θ) * m
    @test [real(mφ.xy); imag(mφ.xy); mφ.z] ≈ vφ
    # Rg
    vg = rotz(φ2) * roty(θ) * rotz(φ1) * v
    mg = Rg(φ1,θ,φ2) * m
    @test [real(mg.xy); imag(mg.xy); mg.z] ≈ vg
    # Rn
    vq = Un(θ, n) * v
    mq = Q(θ, n[1]+n[2]*1im, n[3]) * m
    @test [real(mq.xy); imag(mq.xy); mq.z] ≈ vq

    ## Spinors satify that  |α|^2 + |β|^2 = 1
    @test abs(Rx(θ))                     ≈ [1]
    @test abs(Ry(θ))                     ≈ [1]
    @test abs(Rz(θ))                     ≈ [1]
    @test abs(Rφ(φ,θ))                   ≈ [1]
    @test abs(Q(θ, n[1]+n[2]*1im, n[3])) ≈ [1]

    ## Checking properties of Introduction to the Shinnar-Le Roux algorithm.
    # Rx = Rz(-π/2) * Ry(θ) * Rz(π/2)
    @test rotx(θ) * v ≈ rotz(-π/2) * roty(θ) * rotz(π/2) * v
    @test (Rx(θ) * m).xy ≈ (Rz(-π/2) * Ry(θ) * Rz(π/2) * m).xy
    @test (Rx(θ) * m).z  ≈ (Rz(-π/2) * Ry(θ) * Rz(π/2) * m).z
    # Rφ(φ,θ) = Rz(-φ) Ry(θ) Rz(φ)
    @test (Rφ(φ,θ) * m).xy   ≈ (Rz(-φ) * Ry(θ) * Rz(φ) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ (Rz(-φ) * Ry(θ) * Rz(φ) * m).z
    # Rg(φ1, θ, φ2) = Rz(φ2) Ry(θ) Rz(φ1)
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ (Rz(φ2) * Ry(θ) * Rz(φ1) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ (Rz(φ2) * Ry(θ) * Rz(φ1) * m).z
    # Rg(-φ, θ, φ) = Rz(-φ) Ry(θ) Rz(φ) = Rφ(φ,θ)
    @test rotz(-φ) * roty(θ) * rotz(φ) * v ≈ Un(θ, [sin(φ); cos(φ); 0.0]) * v
    @test (Rg(φ,θ,-φ) * m).xy    ≈ (Rφ(φ,θ) * m).xy
    @test (Rg(φ,θ,-φ) * m).z     ≈ (Rφ(φ,θ) * m).z

    ## Verify trivial identities
    # Rφ is an xy-plane rotation of θ around an axis making an angle of φ with respect to the y-axis
    # Rφ φ=0 = Ry
    @test (Rφ(0,θ) * m).xy   ≈ (Ry(θ) * m).xy
    @test (Rφ(0,θ) * m).z    ≈ (Ry(θ) * m).z
    # Rφ φ=π/2 = Rx
    @test (Rφ(π/2,θ) * m).xy ≈ (Rx(θ) * m).xy
    @test (Rφ(π/2,θ) * m).z  ≈ (Rx(θ) * m).z
    # General rotation Rn
    # Rn n=[1,0,0] = Rx
    @test Un(θ, [1,0,0]) * v ≈ rotx(θ) * v
    @test (Q(θ, 1.0+0.0im, 0.0) * m).xy ≈ (Rx(θ) * m).xy
    @test (Q(θ, 1.0+0.0im, 0.0) * m).z  ≈ (Rx(θ) * m).z
    # Rn n=[0,1,0] = Ry
    @test Un(θ, [0,1,0]) * v ≈ roty(θ) * v
    @test (Q(θ, 0.0+1.0im, 0.0) * m).xy ≈ (Ry(θ) * m).xy
    @test (Q(θ, 0.0+1.0im, 0.0) * m).z  ≈ (Ry(θ) * m).z
    # Rn n=[0,0,1] = Rz
    @test Un(θ, [0,0,1]) * v ≈ rotz(θ) * v
    @test (Q(θ, 0.0+0.0im, 1.0) * m).xy ≈ (Rz(θ) * m).xy
    @test (Q(θ, 0.0+0.0im, 1.0) * m).z  ≈ (Rz(θ) * m).z

    # Associativity
    # Rx
    @test (((Rz(-π/2) * Ry(θ)) * Rz(π/2)) * m).xy ≈ (Rx(θ) * m).xy
    @test (((Rz(-π/2) * Ry(θ)) * Rz(π/2)) * m).z  ≈ (Rx(θ) * m).z
    @test (Rz(-π/2) * (Ry(θ) * (Rz(π/2) * m))).xy ≈ (Rx(θ) * m).xy
    @test (Rz(-π/2) * (Ry(θ) * (Rz(π/2) * m))).z  ≈ (Rx(θ) * m).z
    # Rφ
    @test (Rφ(φ,θ) * m).xy   ≈ (((Rz(-φ) * Ry(θ)) * Rz(φ)) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ (((Rz(-φ) * Ry(θ)) * Rz(φ)) * m).z
    @test (Rφ(φ,θ) * m).xy   ≈ ((Rz(-φ) * (Ry(θ) * Rz(φ))) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ ((Rz(-φ) * (Ry(θ) * Rz(φ))) * m).z
    # Rg
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ (((Rz(φ2) * Ry(θ)) * Rz(φ1)) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ (((Rz(φ2) * Ry(θ)) * Rz(φ1)) * m).z
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ ((Rz(φ2) * (Ry(θ) * Rz(φ1))) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ ((Rz(φ2) * (Ry(θ) * Rz(φ1))) * m).z

    ## Other tests
    # Test Spinor struct
    α, β = rand(2)
    s = Spinor(α, β)
    @test s[1].α ≈ [Complex(α)] && s[1].β ≈ [Complex(β)]
    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", s)
    @test true
end

@testitem "ISMRMRD" tags=[:core, :nomotion] begin
    include("initialize_backend.jl")

    seq = PulseDesigner.EPI_example()[1:10]
    sys = Scanner()
    obj = brain_phantom2D()[1:10]
    parts = kfoldperm(length(obj), 2)

    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "raw"
    sim_params["gpu"] = USE_GPU

    sig1 = simulate(obj[parts[1]], seq, sys; sim_params, verbose=false)
    sig2 = simulate(obj[parts[2]], seq, sys; sim_params, verbose=false)
    sig = simulate(obj, seq, sys; sim_params, verbose=false)

    @test isapprox(sig, sig1 + sig2; rtol=0.001)
end

@testitem "signal_to_raw_data" tags=[:core, :nomotion] begin
    include("initialize_backend.jl")

    seq = PulseDesigner.EPI_example()
    sys = Scanner()
    obj = brain_phantom2D()

    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "mat"
    sim_params["gpu"] = USE_GPU
    sig = simulate(obj, seq, sys; sim_params, verbose=false)

    # Test signal_to_raw_data
    raw = signal_to_raw_data(sig, seq)
    sig_aux = vcat([vec(profile.data) for profile in raw.profiles]...)
    sig_raw = reshape(sig_aux, length(sig_aux), 1)
    @test all(sig .== sig_raw)

    raw = signal_to_raw_data(cat(sig, sig, sig, sig; dims=2), seq)
    head = raw.profiles[1].head
    @test head.available_channels == UInt16(4)
    @test head.active_channels == UInt16(4)
    @test head.channel_mask[1] == UInt64(0xf)
    @test all(iszero, head.channel_mask[2:end])

    seq.DEF["FOV"] = [23e-2, 23e-2, 0]
    raw = signal_to_raw_data(sig, seq)
    sig_aux = vcat([vec(profile.data) for profile in raw.profiles]...)
    sig_raw = reshape(sig_aux, length(sig_aux), 1)
    @test all(sig .== sig_raw)

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", raw)
    @test true
end

@testitem "simulate rejects negative labels" tags=[:core, :nomotion] begin
    obj = Phantom(x=[0.0])
    seq = Sequence([Grad(0, 1e-3)])
    for label in (LabelSet(-1, "LIN"), LabelInc(-1, "LIN"))
        seq.EXT[1] = [label]
        @test_throws Exception simulate(obj, seq, Scanner(); verbose=false)
    end
end

@testitem "simulation precision" tags=[:core, :nomotion] begin
    seq = Sequence()
    @addblock seq += RF([1.0, 2.0, 1.0] .* 1e-6, 1e-4, [0.0, 0.0, 0.0])
    obj = Phantom(
        x=[-5e-4, 0.0, 5e-4], y=[0.0, 4e-4, -4e-4], z=zeros(3),
        ρ=ones(3), T1=fill(1e6, 3), T2=fill(1e6, 3), Δw=zeros(3),
    )

    for (precision, T) in ("f32" => Float32, "f64" => Float64, "bigfloat" => BigFloat)
        sim_params = Dict{String,Any}(
            "gpu" => false,
            "Nthreads" => 1,
            "return_type" => "state",
            "precision" => precision,
            "sim_method" => BlochMagnus1(),
            "Δt" => Inf,
            "Δt_rf" => 50e-6,
        )
        state = simulate(obj, seq, Scanner(); sim_params, verbose=false)
        @test eltype(state.xy) === Complex{T}
        @test eltype(state.z) === T
    end

end

@testitem "Multi-coil simulations" tags=[:core, :nomotion] begin
    coords = [-1.0, 0.0, 1.0] .* 1e-3
    ncoils = 4
    dims = (length(coords), length(coords), length(coords))
    single_coil = fill(1.0 + 0.0im, dims..., 1)
    multi_coil = fill((1.0 + 0.0im) / ncoils, dims..., ncoils)
    single_receiver = ArbitraryCoilSens(coords, coords, coords, single_coil)
    multi_receiver = ArbitraryCoilSens(coords, coords, coords, multi_coil)

    seq = Sequence()
    @addblock seq += RF([1.0, 2.0, 1.0] .* 1e-6, 1e-4, [0.0, 0.0, 0.0])
    @addblock seq += ADC(8, 2e-4)
    obj = Phantom(x=[0.0], y=[0.0], z=[0.0], ρ=[1.0], T1=[1e6], T2=[1e6], Δw=[0.0])

    for sim_method in (BlochSimple(), Bloch())
        @testset "$(nameof(typeof(sim_method)))" begin
            for (precision, T) in ("f32" => Float32, "f64" => Float64)
                @testset "$precision" begin
                    sim_params = Dict{String,Any}(
                        "gpu" => false,
                        "Nthreads" => 1,
                        "return_type" => "mat",
                        "precision" => precision,
                        "sim_method" => sim_method,
                    )
                    raw_uniform = simulate(obj, seq, Scanner(); sim_params, verbose=false)
                    raw_single = simulate(
                        obj, seq, Scanner(receiver=single_receiver); sim_params, verbose=false,
                    )
                    raw_multi = simulate(
                        obj, seq, Scanner(receiver=multi_receiver); sim_params, verbose=false,
                    )

                    @test eltype(raw_uniform) === Complex{T}
                    @test eltype(raw_single) === Complex{T}
                    @test eltype(raw_multi) === Complex{T}
                    @test raw_uniform ≈ raw_single
                    @test raw_uniform ≈ sum(raw_multi; dims=2)
                end
            end
        end
    end

end

@testitem "repeated single-spin FID" tags=[:core, :nomotion] begin
    rf = PulseDesigner.make_block_pulse(π / 2; duration=300e-6, delay=100e-6)
    adc = PulseDesigner.make_adc(4096; dwell=62.5e-6, delay=20e-6)

    seq = Sequence()
    for _ in 1:30
        @addblock seq += (rf, Duration(20e-3))
        @addblock seq += (adc, Duration(10.0))
    end

    obj = Phantom(x=[0.0], y=[0.0], z=[0.0], T1=[1.0], T2=[0.1], Δw=[100.0])
    sim_params = Dict{String,Any}(
        "gpu" => false,
        "return_type" => "mat",
        "sim_method" => Bloch(),
    )
    raw = simulate(obj, seq, Scanner(); sim_params, verbose=false)
    first_adc = raw[1:4096:end, 1, 1]

    # Regression: timing roundoff used to place RF samples incorrectly,
    # changing the effective flip angle and the acquired signal in later TRs.
    @test minimum(abs, first_adc) / maximum(abs, first_adc) > 0.99
    @test maximum(abs, first_adc ./ first(first_adc) .- 1) < 1e-3
end

@testitem "Bloch" tags=[:important, :core, :nomotion, :bloch] begin
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_sphere_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    for sim_method in (
        KomaMRICore.Bloch(),
        KomaMRICore.BlochMagnusConst1(),
        KomaMRICore.BlochMagnusLin2(),
        KomaMRICore.BlochMagnusLinComm2(),
        KomaMRICore.BlochMagnusQuad2(),
        KomaMRICore.BlochMagnusQuad4(),
        KomaMRICore.BlochMagnusGL2(),
        KomaMRICore.BlochMagnusGL4()
    )
        @testset "$(nameof(typeof(sim_method)))" begin
            sim_params = Dict{String, Any}(
                "gpu"=>USE_GPU,
                "sim_method"=>sim_method,
                "return_type"=>"mat"
            )
            sig = simulate(obj, seq, sys; sim_params, verbose=false)
            sig = sig / prod(size(obj))

            @test NRMSE(sig, sig_jemris) < 1 #NRMSE < 1%
        end
    end

    if USE_GPU
        coords = [-1.0, 0.0, 1.0] .* 1e-3
        ncoils = 4
        coil_sens = fill(ComplexF32(1 / ncoils), length(coords), length(coords), length(coords), ncoils)
        receivers = (
            BirdcageCoilSens(; ncoils),
            ArbitraryCoilSens(coords, coords, coords, coil_sens),
        )
        multi_coil_obj = Phantom(
            x=[-5e-4, 0.0, 5e-4], y=[0.0, 4e-4, -4e-4], z=zeros(3),
            ρ=ones(3), T1=fill(1e6, 3), T2=fill(1e6, 3), Δw=zeros(3),
        )
        multi_coil_seq = Sequence()
        @addblock multi_coil_seq += RF([1.0, 2.0, 1.0] .* 1e-6, 1e-4) + ADC(4, 1e-4)
        @addblock multi_coil_seq += ADC(4, 2e-4)
        multi_coil_methods = (
            Bloch(), BlochMagnusConst1(), BlochMagnusLin2(), BlochMagnusMid2(),
            BlochMagnusLinComm2(), BlochMagnusQuad2(), BlochMagnusQuad4(),
            BlochMagnusGL2(), BlochMagnusGL4(), BlochMagnusBGL4(), BlochMagnusBGL6(),
        )

        for receiver in receivers, sim_method in multi_coil_methods
            sim_params = Dict{String,Any}(
                "gpu" => false,
                "Nthreads" => 1,
                "precision" => "f32",
                "return_type" => "mat",
                "sim_method" => sim_method,
            )
            reference = simulate(
                multi_coil_obj, multi_coil_seq, Scanner(; receiver);
                sim_params, verbose=false,
            )
            sim_params["gpu"] = true
            signal = simulate(
                multi_coil_obj, multi_coil_seq, Scanner(; receiver);
                sim_params, verbose=false,
            )
            @test signal ≈ reference
        end
    end
end

@testitem "Bloch waveform event type accuracy" tags=[:core, :nomotion] begin
    using OrdinaryDiffEqTsit5
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    Tpulse = 1e-3
    Tgrad = 1e-3
    Tadc = 1e-3
    Nadc = 6
    M0 = 1.0
    T1 = 1000e-3
    T2 = 40e-3
    Δw = 2π * 30
    x0 = 1e-2
    B1 = 1.5e-6 * cis(π / 7)
    Gx = 0.2e-3

    sys = Scanner()
    obj = Phantom(x=[x0], ρ=[M0], T1=[T1], T2=[T2], Δw=[Δw])

    grad_events = (
        "trap" => Grad(Gx, Tgrad),
        "uniform" => Grad([Gx, Gx], Tgrad),
        "time-shaped" => Grad([Gx, Gx], [Tgrad]),
    )
    rf_events = (
        "block" => RF(B1, Tpulse),
        "uniform" => RF([B1, B1], Tpulse),
        "time-shaped" => RF([B1, B1], [Tpulse]),
    )

    function waveform_sequence(grad, rf)
        seq = Sequence()
        @addblock seq += rf + (ADC(Nadc, Tadc), x=grad)
        return seq
    end

    ref_seq = waveform_sequence(grad_events[1][2], rf_events[1][2])
    ref_adc_times = get_adc_sampling_times(ref_seq)
    mxy_diffeq = diffeq_signal(ref_seq, obj; tstops=[Tpulse, Tpulse + Tgrad])

    for (grad_name, grad) in grad_events, (rf_name, rf) in rf_events
        seq = waveform_sequence(grad, rf)
        @test get_adc_sampling_times(seq) ≈ ref_adc_times
        for sim_method in (Bloch(), BlochMagnusConst1(), BlochMagnusLin2(), BlochMagnusLinComm2(), BlochMagnusQuad2(), BlochMagnusQuad4(), BlochMagnusGL2(), BlochMagnusGL4())
            @testset "$grad_name/$rf_name $(nameof(typeof(sim_method)))" begin
                sim_params = Dict{String, Any}(
                    "gpu" => USE_GPU,
                    "return_type" => "mat",
                    "sim_method" => sim_method,
                )
                raw = simulate(obj, seq, sys; sim_params, verbose=false)[:, 1, 1]
                @test NRMSE(raw, mxy_diffeq) < 0.1
            end
        end
    end
end

@testitem "Bloch_RF_accuracy" tags=[:important, :core, :nomotion] begin
    using OrdinaryDiffEqTsit5
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    # Seq params
    Tadc = 1e-3
    Trf = Tadc
    B1 = 2e-6 * (Tadc / Trf)
    rf_phase = [0, π/2]
    Gx = 25e-3
    Nadc = 6
    # Phantom params
    M0 = 1.0
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    x = 1e-2
    
    ## Solving using KomaMRI
    seq = Sequence()
    @addblock seq += ADC(Nadc, Tadc)
    @addblock seq += (RF(B1 .* cis(rf_phase[1]), Trf), x=Grad(Gx, Trf))
    @addblock seq += ADC(Nadc, Tadc)
    # This introduces an RF-ADC overlap!!!
    @addblock seq += (RF(B1 .* cis(rf_phase[2]), Trf), ADC(Nadc, 2Tadc, Trf/2))

    sys = Scanner()
    obj = Phantom(x = [x], ρ = [M0], T1 = [T1], T2 = [T2], Δw = [Δw])

    mxy_diffeq = diffeq_signal(seq, obj)

    ## Solve with KomaMRI
    methods_to_test = [Bloch(), BlochMagnusConst1(), BlochMagnusLin2(), BlochMagnusLinComm2(), BlochMagnusQuad2(), BlochMagnusQuad4(), BlochMagnusGL2(), BlochMagnusGL4()]
    sim_params_to_test = [Dict{String, Any}("Δt_rf"=>1e-5, "return_type"=>"mat", "sim_method"=>method) for method in methods_to_test]
    for sim_params in sim_params_to_test
        @testset "$(sim_params["sim_method"])" begin 
            raw = simulate(obj, seq, sys; sim_params, verbose=false)[:, 1, 1]
            @test NRMSE(raw, mxy_diffeq) < 0.1
        end
    end
end

@testitem "Bloch_phase_compensation" tags=[:important, :core, :nomotion] begin
    include("initialize_backend.jl")

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom(x=[0.],T1=[T1],T2=[T2],Δw=[Δw])

    rf_phase = 2π*rand()
    seq1 = Sequence()
    seq1 += RF(B1, Trf)
    seq1 += ADC(N, Tadc)

    seq2 = Sequence()
    seq2 += RF(B1 .* exp(1im*rf_phase), Trf)
    seq2 += ADC(N, Tadc, 0, 0, rf_phase)

    for sim_method in (
        KomaMRICore.Bloch(),
        KomaMRICore.BlochMagnusConst1(),
        KomaMRICore.BlochMagnusLin2(),
        KomaMRICore.BlochMagnusLinComm2(),
        KomaMRICore.BlochMagnusQuad2(),
        KomaMRICore.BlochMagnusQuad4(),
        KomaMRICore.BlochMagnusGL2(),
        KomaMRICore.BlochMagnusGL4()
    )
        @testset "$(nameof(typeof(sim_method)))" begin
            sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>USE_GPU, "sim_method"=>sim_method)
            raw1 = simulate(obj, seq1, sys; sim_params, verbose=false)
            raw2 = simulate(obj, seq2, sys; sim_params, verbose=false)

            @test raw1.profiles[1].data ≈ raw2.profiles[1].data
        end
    end
end

@testitem "BlochDict" tags=[:important, :core, :nomotion, :blochdict] begin
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    seq = seq_epi_100x100_TE100_FOV230()
    obj = Phantom(x=[0.], T1=[1000e-3], T2=[100e-3])
    sys = Scanner()
    sim_params = Dict(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat")
    sig = simulate(obj, seq, sys; sim_params, verbose=false)
    sig = sig / prod(size(obj))
    sim_params["sim_method"] = KomaMRICore.BlochDict()
    sig2 = simulate(obj, seq, sys; sim_params, verbose=false)
    sig2 = sig2 / prod(size(obj))
    @test sig ≈ sig2

    sig_jemris = signal_sphere_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()
    sig = simulate(obj, seq, sys; sim_params, verbose=false)
    sig = sum(sig; dims=2) / prod(size(obj))
    @test NRMSE(sig, sig_jemris) < 1 #NRMSE < 1%

    sim_params["sim_method"] = KomaMRICore.BlochDict(save_Mz=true)
    sig2 = simulate(obj[1], seq[1:100], sys; sim_params, verbose=false)
    @test true # Just checking that it runs, TODO: compare to DiffEq

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", KomaMRICore.BlochDict())
    @test true
end

@testitem "BlochSimple" tags=[:important, :core, :nomotion, :blochsimple] begin
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_sphere_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = simulate(obj, seq, sys; sim_params, verbose=false)
    sig = sig / prod(size(obj))

    @test NRMSE(sig, sig_jemris) < 1 #NRMSE < 1%
end

@testitem "simulate_slice_profile" tags=[:core, :nomotion] begin
    include("initialize_backend.jl")

    # This is a sequence with a sinc RF 30° excitation pulse
    sys = Scanner()
    sys.limits.Smax = 50
    B1 = 4.92e-6
    Trf = 3.2e-3
    zmax = 2e-2
    fmax = 5e3
    z = range(-zmax, zmax, 400)
    Gz = fmax / (γ * zmax)
    seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0; 0; Gz], TBP=8)

    sim_params = Dict{String, Any}(
        "Δt_rf" => Trf / length(seq.RF.A[1]),
        "gpu" => USE_GPU)

    @testset "frequency offset shifts sinc slice profile" begin
        sample_shift = 40
        Δf = γ * Gz * sample_shift * step(z)
        shifted_seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0; 0; Gz], Δf, TBP=8)
        methods_to_test = (
            KomaMRICore.Bloch(),
            KomaMRICore.BlochDict(),
            KomaMRICore.BlochMagnusConst1(),
            KomaMRICore.BlochMagnusLin2(),
            KomaMRICore.BlochMagnusLinComm2(),
            KomaMRICore.BlochMagnusQuad2(),
            KomaMRICore.BlochMagnusQuad4(),
            KomaMRICore.BlochMagnusGL2(),
            KomaMRICore.BlochMagnusGL4(),
        )
        for sim_method in methods_to_test
            @testset "$(nameof(typeof(sim_method)))" begin
                for max_rf_block_length in (Inf, 30, 1)
                    @testset "max_rf_block_length=$max_rf_block_length" begin
                        shifted_sim_params = copy(sim_params)
                        shifted_sim_params["sim_method"] = sim_method
                        shifted_sim_params["max_rf_block_length"] = max_rf_block_length
                        base_sim_params = copy(shifted_sim_params)
                        M_base = simulate_slice_profile(seq; z, sim_params=base_sim_params, verbose=false)
                        M_shifted = simulate_slice_profile(shifted_seq; z, sim_params=shifted_sim_params, verbose=false)

                        profile = abs.(M_base.xy)
                        shifted_profile = abs.(M_shifted.xy)
                        expected = profile[1:(end - sample_shift)]
                        shifted = shifted_profile[(1 + sample_shift):end]
                        @test shifted ≈ expected
                    end
                end
            end
        end
    end
end

@testitem "freq_in_phase preserves RF frequency offsets" tags=[:core, :nomotion, :gpu] begin
    include("initialize_backend.jl")

    sim_params = Dict{String,Any}(
        "gpu" => USE_GPU,
        "return_type" => "state",
        "Δt" => Inf,
        "Δt_rf" => 5e-6,
    )
    methods = (Bloch(), BlochSimple(), BlochMagnus1(), BlochMagnus2(), BlochMagnus4(), BlochMagnus6())

    @testset "slice-selective excitation" begin
        sys = Scanner()
        sys.limits.Smax = 50
        B1 = 4.92e-6
        Trf = 3.2e-3
        z0 = 4e-3
        z = [-z0, 0.0, z0]
        Gz = 5e3 / (γ * 0.02)
        Δf = γ * Gz * z0
        seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0; 0; Gz], TBP=8)
        shifted_seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0; 0; Gz], Δf=Δf, TBP=8)

        for method in methods, freq_in_phase in (false, true)
            params = copy(sim_params)
            params["sim_method"] = method
            params["freq_in_phase"] = freq_in_phase
            profile = abs.(simulate_slice_profile(seq; z, sim_params=copy(params), verbose=false).xy)
            shifted_profile = abs.(simulate_slice_profile(shifted_seq; z, sim_params=copy(params), verbose=false).xy)
            @test z[argmax(profile)] == 0.0
            @test z[argmax(shifted_profile)] == z0
        end
    end

    @testset "adiabatic inversion" begin
        duration = 18.3e-3
        b1max = 13.5e-6
        β̂ = 4
        μ = 6
        β = 2β̂ / duration
        t = range(-duration / 2, duration / 2, 4001)
        B1 = b1max .* sech.(β .* t)
        Δf = -μ * β .* tanh.(β .* t) ./ (2π)
        shift = 1000.0
        freqs = [-shift, 0.0, shift]
        obj = Phantom(x=zeros(length(freqs)), y=zeros(length(freqs)), z=zeros(length(freqs)), Δw=2π .* freqs)
        seq = Sequence()
        shifted_seq = Sequence()
        @addblock seq += RF(B1, duration, Δf, 0)
        @addblock shifted_seq += RF(B1, duration, Δf .+ shift, 0)

        for method in methods, freq_in_phase in (false, true)
            params = copy(sim_params)
            params["sim_method"] = method
            params["freq_in_phase"] = freq_in_phase
            profile = simulate(obj, seq, Scanner(); sim_params=copy(params), verbose=false).z
            shifted_profile = simulate(obj, shifted_seq, Scanner(); sim_params=copy(params), verbose=false).z
            @test freqs[argmin(profile)] == 0.0
            @test freqs[argmin(shifted_profile)] == shift
        end
    end

    @testset "fat saturation" begin
        Δf = -440.0
        seq = Sequence()
        @addblock seq += RF(fill(1.0e-6, 65), 1e-3, Δf)
        seq.RF[1].A .*= 90 / only(get_flip_angles(seq))
        freqs = [-Δf, 0.0, Δf]
        obj = Phantom(x=zeros(length(freqs)), Δw=2π .* freqs)

        for method in methods, freq_in_phase in (false, true)
            params = copy(sim_params)
            params["sim_method"] = method
            params["freq_in_phase"] = freq_in_phase
            state = simulate(obj, seq, Scanner(); sim_params=params, verbose=false)
            @test freqs[argmin(abs.(state.z))] == Δf
        end
    end
end

@testitem "GPU Functions" tags=[:core, :nomotion, :gpu] begin
    using Suppressor
    import KernelAbstractions as KA
    include("initialize_backend.jl")

    x = ones(Float32, 1000)

    begin
        if USE_GPU
            y = x |> gpu
            @test KA.get_backend(y) isa KA.GPU
            y = y |> cpu
            @test KA.get_backend(y) isa KA.CPU
        else
            # Test that gpu and cpu are no-ops
            y = x |> gpu
            @test y == x
            y = y |> cpu
            @test y == x
        end
    end

    @suppress print_devices()
    @test true
end

# --------- Motion-related tests -------------
# We compare with the result given by OrdinaryDiffEqTsit5
@testitem "Motion" tags=[:core, :motion] begin
    using OrdinaryDiffEqTsit5
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    Nadc = 25
    M0 = 1.0
    T1 = 100e-3
    T2 = 10e-3
    B1 = 20e-6
    Trf = 3e-3
    φ = π / 4
    duration = 2*Trf

    Gx = 1e-3
    Gy = 1e-3
    Gz = 0.0

    motions = [
        translate(0.1, 0.1, 0.0, TimeRange(0.0, 1.0)),
        rotate(0.0, 0.0, 45.0, TimeRange(0.0, 1.0)),
        heartbeat(-0.6, 0.0, 0.0, Periodic(period=1.0)),
        path([0.0 0.0], [0.0 1.0], [0.0 0.0], TimeRange(0.0, 10.0)),
        flowpath([0.0 0.0], [0.0 1.0], [0.0 0.0], [false false], TimeRange(0.0, 10.0)) # We should find a way to test this when spin_reset flags are true
    ]

    x0 = [0.1]
    y0 = [0.1]
    z0 = [0.0]

    for sim_method in [
        KomaMRICore.Bloch(),
        KomaMRICore.BlochSimple(),
        KomaMRICore.BlochDict(),
        KomaMRICore.BlochMagnusConst1(),
        KomaMRICore.BlochMagnusLin2(),
        KomaMRICore.BlochMagnusLinComm2(),
        KomaMRICore.BlochMagnusQuad2(),
        KomaMRICore.BlochMagnusQuad4(),
        KomaMRICore.BlochMagnusGL2(),
        KomaMRICore.BlochMagnusGL4()
    ]
        @testset "$(typeof(sim_method))" begin
            for motion in motions
                seq = Sequence()
                @addblock seq += RF(cis(φ) .* B1, Trf) + (
                    ADC(Nadc, duration - Trf, Trf),
                    x=Grad(Gx, duration),
                    y=Grad(Gy, duration),
                    z=Grad(Gz, duration),
                )
                obj = Phantom(x = x0, y = y0, z = z0, ρ = [M0], T1 = [T1], T2 = [T2], motion = motion)
                sys = Scanner()
                mxy_diffeq = diffeq_signal(seq, obj)
                sim_params = Dict{String, Any}(
                    "sim_method"=>sim_method,
                    "return_type"=>"mat",
                    "gpu" => USE_GPU
                )
                raw_aux = simulate(obj, seq, sys; sim_params, verbose=false)
                raw = raw_aux[:, 1, 1]
                @test NRMSE(raw, mxy_diffeq) < 1
            end
        end
    end
end
