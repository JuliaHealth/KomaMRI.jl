using FiniteDifferences: central_fdm, grad

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
        "Δt" => 50e-6,
        "Δt_rf" => 25e-6,
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
    seq_aux = copy(params.seq)
    seq_aux.RF[1].A = complex.(rf_samples) .* params.rf_scale
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
    seq = Sequence(
        params.seq.GR,
        Reactant.to_rarray.(params.seq.RF),
        params.seq.ADC,
        params.seq.DUR,
        params.seq.EXT,
        params.seq.DEF,
    )
    return merge(params, (;
        seq,
        obj=Reactant.to_rarray(params.obj),
        target_profile=Reactant.to_rarray(params.target_profile),
    ))
end
