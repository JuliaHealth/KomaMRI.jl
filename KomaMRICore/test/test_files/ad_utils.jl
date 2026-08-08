using FiniteDifferences: central_fdm, grad

const BLOCHSIMPLE_PARALLEL_AD_RF0 = [1.3, 1.7, 1.1]

function blochsimple_parallel_ad_vector(template, values, ::Type{T}=Float64) where {T}
    out = similar(template, T, length(values))
    Reactant.allowscalar() do
        for i in eachindex(values)
            out[i] = values[i]
        end
    end
    return out
end

function blochsimple_parallel_ad_sequence(rf, t_start, t_end; excitation)
    n = length(rf)
    t = collect(range(t_start, t_end; length=n))
    z = zero.(rf)
    B1 = excitation ? complex.(rf) .* 1e-6 : complex.(z)
    return DiscreteSequence(
        z,
        copy(z),
        copy(z),
        B1,
        copy(z),
        copy(z),
        falses(n),
        fill(excitation, n - 1),
        t,
        diff(t),
    )
end

function blochsimple_parallel_ad_loss(rf)
    obj = Phantom(
        x=blochsimple_parallel_ad_vector(rf, [0.0, 1e-2]),
        ρ=blochsimple_parallel_ad_vector(rf, [1.0, 0.75]),
        T1=blochsimple_parallel_ad_vector(rf, [1.0, 0.8]),
        T2=blochsimple_parallel_ad_vector(rf, [0.08, 0.12]),
        Δw=blochsimple_parallel_ad_vector(rf, 2π .* [10.0, -12.0]),
    )
    M = Mag(
        blochsimple_parallel_ad_vector(rf, zeros(2), ComplexF64),
        blochsimple_parallel_ad_vector(rf, [1.0, 0.75]),
    )
    signal = similar(rf, ComplexF64, 0, 1, 1)
    method = BlochSimple()
    backend = KomaMRICore.KA.CPU()
    prealloc = KomaMRICore.DefaultPrealloc()

    KomaMRICore.run_spin_excitation_parallel!(
        obj,
        blochsimple_parallel_ad_sequence(rf, 0.0, 0.6e-3; excitation=true),
        signal,
        M,
        method,
        256,
        backend,
        prealloc;
        Nthreads=1,
    )
    KomaMRICore.run_spin_precession_parallel!(
        obj,
        blochsimple_parallel_ad_sequence(rf, 0.6e-3, 1.4e-3; excitation=false),
        signal,
        M,
        method,
        256,
        backend,
        prealloc;
        Nthreads=1,
    )

    target = blochsimple_parallel_ad_vector(rf, [0.4, 0.2])
    return sum(abs2, M.xy) + sum(abs2, M.z .- target)
end

blochsimple_parallel_ad_fd_gradient(rf=BLOCHSIMPLE_PARALLEL_AD_RF0) =
    grad(central_fdm(5, 1), blochsimple_parallel_ad_loss, rf)[1]

const BLOCH_NODE_AD_X0 = [0.0, 3.5, 0.0, 1.0, 0.75, 0.5]

function bloch_node_ad_parameters(sim_method)
    rf_duration = 0.6e-3
    total_duration = 1.4e-3
    z = collect(range(-4e-3, 4e-3; length=3))
    obj = Phantom(
        x=zeros(length(z)),
        y=zeros(length(z)),
        z=z,
        ρ=ones(length(z)),
        T1=ones(length(z)),
        T2=fill(0.1, length(z)),
        Δw=zeros(length(z)),
    )
    rf = RF(
        zeros(ComplexF64, 7),
        rf_duration,
        0.0,
        0.0,
        rf_duration / 2,
        0.0,
        Excitation(),
        Val(:preserve),
    )
    zero_gradient = Grad(0.0, 0.0)
    slice_gradient = Grad(8e-3, rf_duration)
    seq = Sequence(
        reshape([zero_gradient, zero_gradient, slice_gradient], 3, 1),
        reshape([rf], 1, 1),
        [ADC(3, total_duration, 0.0)],
        [total_duration],
        [Extension[]],
        Dict{String,Any}(),
    )
    sim_params = Dict{String,Any}(
        "sim_method" => sim_method,
        "gpu" => false,
        "Nthreads" => 1,
        "return_type" => "mat",
        "precision" => "f64",
        "sampling_rule" => MaxStepSizeRule(50e-6, 25e-6),
    )
    params = (;
        seq,
        obj,
        sys=Scanner(),
        sim_params,
        target_profile=zeros(ComplexF64, 3),
        node_times=range(0.0, rf_duration; length=3),
        rf_times=range(0.0, rf_duration; length=7),
        rf_scale=1e-6,
    )
    target_profile = copy(bloch_node_ad_forward([0.0, 5.0, 0.0, 0.9, 0.7, 0.6], params))
    return merge(params, (; target_profile))
end

function bloch_node_ad_forward(x, params)
    rf_samples = KomaMRIBase.linear_interpolate_samples(
        (t=params.node_times, A=x[1:3]),
        params.rf_times,
    )
    seq_aux = KomaMRIBase.set_rf_amplitude(
        params.seq,
        complex.(rf_samples) .* params.rf_scale,
    )
    obj_aux = copy(params.obj)
    obj_aux.ρ .= x[4:6]
    return simulate(
        obj_aux,
        seq_aux,
        params.sys;
        sim_params=params.sim_params,
        verbose=false,
    )
end

function bloch_node_ad_loss(x, params)
    signal = bloch_node_ad_forward(x, params)
    return sum(abs2, signal .- params.target_profile) / length(signal)
end

bloch_node_ad_fd_gradient(params, x=BLOCH_NODE_AD_X0) =
    grad(central_fdm(5, 1), x -> bloch_node_ad_loss(x, params), x)[1]

function bloch_node_ad_reactant_parameters(params)
    return merge(params, (;
        obj=Reactant.to_rarray(params.obj),
        target_profile=Reactant.to_rarray(params.target_profile),
    ))
end
