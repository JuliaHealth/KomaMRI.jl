@testitem "Scanner" tags=[:files] begin
    using KomaMRIBase

    limits = HardwareLimits(
        B0=3.0, B1=20e-6, Gmax=80e-3, Smax=200.0,
        ADC_Δt=1e-6, DUR_Δt=2e-5, GR_Δt=2e-5, RF_Δt=2e-6,
        RF_ring_down_time=30e-6, RF_dead_time=120e-6, ADC_dead_time=15e-6,
    )
    x, y = [-0.08, 0.02, 0.11], [-0.04, 0.07]
    planar = ArbitraryCoilSens(x, y, [0.0],
        ComplexF32[complex(i + j + k, c - i) for i in 1:3, j in 1:2, k in 1:1, c in 1:2])
    volume = ArbitraryCoilSens(x, y, [-0.03, 0.05],
        ComplexF64[complex(i + j + k, c - i) for i in 1:3, j in 1:2, k in 1:2, c in 1:3])

    mktempdir() do dir
        for receiver in (UniformCoilSens(), BirdcageCoilSens(ncoils=4, radius=0.18, L=0.24), planar, volume)
            original = Scanner(; limits, receiver)
            filename = joinpath(dir, "scanner.sys")
            write_scanner(original, filename)
            restored = read_scanner(filename)

            # The selected models and every hardware/model parameter survive the file round trip.
            for name in fieldnames(Scanner)
                expected, actual = getproperty(original, name), getproperty(restored, name)
                @test nameof(typeof(actual)) == nameof(typeof(expected))
                @test all(getproperty(actual, field) == getproperty(expected, field)
                    for field in fieldnames(typeof(expected)))
            end

            # Restoring a map preserves complex coil responses at grid, interpolated, and exterior positions.
            positions = ([-0.08, 0.0, 0.5], [-0.04, 0.01, 0.5], [0.0, 0.0, 0.5])
            @test get_sens(restored.receiver, positions...) == get_sens(receiver, positions...)
        end
    end
end
