@testitem "Gradient spectrum MATLAB Pulseq parity" tags=[:koma] begin
    using MAT

    fixtures = joinpath(@__DIR__, "test_files", "gradient_spectrum")
    reference = matread(joinpath(fixtures, "pulseq_reference.mat"))
    for name in ("spiral", "epi_ramp")
        seq = read_seq(joinpath(pkgdir(KomaMRI), "examples", "1.sequences", name * ".seq"); verbose=false)
        for fmax in (3000, 1200)
            @testset "$name, $fmax Hz" begin
                expected = reference[name]["f$fmax"]
                spectrum = fmax == 3000 ? KomaMRI.gradient_spectrum(seq) : KomaMRI.gradient_spectrum(seq; fmax)
                # Float64 interpolation and the 15,000-point FFT accumulate O(eps * log2(N)) roundoff.
                roundoff = 128eps(Float64) * log2(15000)

                # MATLAB's unnormalised axis and RSS spectra agree for short and overlapping windows.
                @test spectrum.F ≈ vec(expected["F"]) rtol=roundoff
                for axis in axes(spectrum.Rax, 1)
                    @test spectrum.Rax[axis, :] ≈ expected["Rax"][axis, :] rtol=roundoff
                end
                @test spectrum.R ≈ vec(expected["R"]) rtol=roundoff
            end
        end
    end
end
