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
