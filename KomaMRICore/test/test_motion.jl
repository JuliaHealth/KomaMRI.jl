@testitem "Bloch SimpleAction" begin 
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_simple_motion()
    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE SimpleAction: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch ArbitraryAction" begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_arbitrary_motion()
    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE ArbitraryAction: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "BlochSimple SimpleAction" begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_simple_motion()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE SimpleAction BlochSimple: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "BlochSimple ArbitraryAction" begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_arbitrary_motion()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE ArbitraryAction BlochSimple: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end