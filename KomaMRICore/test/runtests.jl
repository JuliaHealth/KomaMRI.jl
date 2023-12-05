using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:core in ti.tags) #verbose=true

# Test ISMRMRD
@testitem "ISMRMRD" begin
    using JLD2

    path = joinpath(@__DIR__, "test_files")
    seq = load_object(joinpath(path, "radial_JEMRIS.jld2"))

    # Test ISMRMRD
    raw = load_object(joinpath(path, "Koma_signal.jld2"))
    @test raw.params["protocolName"] == "epi"
    @test raw.params["institutionName"] == "Pontificia Universidad Catolica de Chile"
    @test raw.params["encodedSize"] ≈ [101, 101, 1]
    @test raw.params["reconSize"] ≈ [102, 102, 1]
    @test raw.params["patientName"] == "brain2D_axial"
    @test raw.params["trajectory"] == "other"
    @test raw.params["systemVendor"] == "KomaMRI.jl"

    # Test signal_to_raw_data
    sig_aux = vcat([vec(profile.data) for profile in raw.profiles]...)
    sig = reshape(sig_aux, length(sig_aux), 1)
    rawmrd = signal_to_raw_data(sig, seq)
    @test rawmrd.params["institutionName"] == raw.params["institutionName"]

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", rawmrd)
    @test true
end

@testitem "Bloch_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor, JLD2, KomaMRIBase

    path = joinpath(@__DIR__, "test_files")
    seq = load_object(joinpath(path, "epi_100x100_TE100_FOV230.jld2"))
    obj = Phantom(; load_object(joinpath(path, "sphere_chemical_shift.jld2"))...)
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>false,
        "Nthreads"=>1,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    sig_jemris = load_object(joinpath(path, "jemris_signals_epi_sphere_cs.jld2"))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_CPU_multi_thread" tags=[:important, :core] begin
    using Suppressor, JLD2, KomaMRIBase

    path = joinpath(@__DIR__, "test_files")
    seq = load_object(joinpath(path, "epi_100x100_TE100_FOV230.jld2"))
    obj = Phantom(; load_object(joinpath(path, "sphere_chemical_shift.jld2"))...)
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>false,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    sig_jemris = load_object(joinpath(path, "jemris_signals_epi_sphere_cs.jld2"))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_GPU" tags=[:important, :skipci, :core] begin
    using Suppressor, JLD2, KomaMRIBase

    path = joinpath(@__DIR__, "test_files")
    seq = load_object(joinpath(path, "epi_100x100_TE100_FOV230.jld2"))
    obj = Phantom(; load_object(joinpath(path, "sphere_chemical_shift.jld2"))...)
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>true,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    sig_jemris = load_object(joinpath(path, "jemris_signals_epi_sphere_cs.jld2"))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_CPU_RF_accuracy_single_thread" tags=[:important, :core] begin
    using Suppressor, KomaMRIBase

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false, "Nthreads"=>1)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "Bloch_CPU_RF_accuracy_multi_thread" tags=[:important, :core] begin
    using Suppressor, KomaMRIBase

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "Bloch_GPU_RF_accuracy" tags=[:important, :core] begin
    using Suppressor, KomaMRIBase

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>true)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "BlochDict_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor, JLD2, KomaMRIBase

    path = joinpath(@__DIR__, "test_files")
    seq = load_object(joinpath(path, "epi_100x100_TE100_FOV230.jld2"))
    obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
    sys = Scanner()
    sim_params = Dict("gpu"=>false, "Nthreads"=>1, "sim_method"=>KomaMRICore.Bloch(), "return_type"=>"mat")
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    sim_params["sim_method"] = KomaMRICore.BlochDict()
    sig2 = simulate(obj, seq, sys; sim_params)
    sig2 = sig2 / prod(size(obj))
    @test sig ≈ sig2

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", KomaMRICore.BlochDict())
    @test true
end
