using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:koma in ti.tags) #verbose=true

# include("../KomaMRICore/test/runtests.jl")
# include("../KomaMRIPlots/test/runtests.jl")

@testitem "MRIReco recon" tags=[:koma] begin
    #Sanity check 1
    A = rand(5,5,3)
    B = KomaMRI.fftc(KomaMRI.ifftc(A))
    @test A ≈ B
    
    #Sanity check 2
    B = KomaMRI.ifftc(KomaMRI.fftc(A))
    @test A ≈ B

    #MRIReco.jl
    path = @__DIR__
    fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
    raw = RawAcquisitionData(fraw)
    acq = AcquisitionData(raw)

    @testset "MRIReco_direct" begin
        Nx, Ny = raw.params["reconSize"][1:2]
        recParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx,Ny), :densityWeighting=>true)
        img = reconstruction(acq, recParams)
        @test true                #If the previous line fails the test will fail
    end

    #Test MRIReco regularized recon (with a λ)
    @testset "MRIReco_standard" begin
        #???
    end

end

@testitem "KomaCLI" tags=[:koma] begin
    using KomaMRI

    @testset "Defaults" begin
        opts = KomaMRI.parse_cli_args(String[])
        @test isnothing(opts.sequence)
        @test isnothing(opts.phantom)
        @test isnothing(opts.scanner)
        @test isnothing(opts.sim_output)
        @test isnothing(opts.recon_output)
        @test opts.backend == "CPU"
        @test opts.sim_params["gpu"] == false
        KomaMRI.load_cli_backend!(opts)
        @test opts.sim_params["gpu"] == false
        @test opts.recon_params[:reco] == "direct"
    end

    @testset "Inputs and outputs" begin
        opts = KomaMRI.parse_cli_args(["-i", "obj.phantom", "seq.seq", "scanner.sys", "-o", "raw.mrd", "image.mat"])
        @test opts.sequence == "seq.seq"
        @test opts.phantom == "obj.phantom"
        @test opts.scanner == "scanner.sys"
        @test opts.sim_output == "raw.mrd"
        @test opts.recon_output == "image.mat"

        opts = KomaMRI.parse_cli_args(["--inputs", "seq.seq", "obj.h5", "--outputs", "_", "image.mat"])
        @test opts.sequence == "seq.seq"
        @test opts.phantom == "obj.h5"
        @test isnothing(opts.sim_output)
        @test opts.recon_output == "image.mat"
    end

    @testset "Parameters and backend" begin
        opts = KomaMRI.parse_cli_args(["--backend=CPU", "-s", "gpu=true", "-s", "Nthreads=4", "-s", "sim_method=BlochMagnus4", "-r", "reco=direct", "-r", "shape=(2,3)"])
        @test opts.backend == "CPU"
        @test opts.sim_params["gpu"] == true
        @test opts.sim_params["Nthreads"] == 4
        @test opts.sim_params["sim_method"] isa KomaMRI.BlochMagnus4
        @test opts.recon_params[:reco] == "direct"
        @test opts.recon_params[:shape] == (2, 3)

        KomaMRI.load_cli_backend!(opts)
        @test opts.sim_params["gpu"] == false
    end

    @testset "Preferences" begin
        opts = KomaMRI.merge_cli_preferences!(
            KomaMRI.CLIOptions(),
            Dict{String,Any}(
                "backend" => "Metal",
                "inputs" => Dict{String,Any}("sequence" => "mysequence.seq", "phantom" => "myphantom.phantom", "scanner" => "scanner.sys"),
                "outputs" => Dict{String,Any}("rawdata" => "raw.mrd", "image" => "image.mat"),
                "sim_params" => Dict{String,Any}("sim_method" => "BlochMagnus4", "precision" => "f32"),
                "recon_params" => Dict{String,Any}("reco" => "direct"),
            ),
        )
        @test opts.sequence == "mysequence.seq"
        @test opts.phantom == "myphantom.phantom"
        @test opts.scanner == "scanner.sys"
        @test opts.sim_output == "raw.mrd"
        @test opts.recon_output == "image.mat"
        @test opts.backend == "Metal"
        @test opts.sim_params["sim_method"] isa KomaMRI.BlochMagnus4
        @test opts.sim_params["precision"] == "f32"
        @test opts.recon_params[:reco] == "direct"

        opts = KomaMRI.parse_cli_args(["-b", "CPU", "-s", "precision=f64"], opts)
        @test opts.backend == "CPU"
        @test opts.sim_params["precision"] == "f64"
    end

    @testset "Errors" begin
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs", "a.seq", "b.phantom", "c.sys", "d.seq"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs", "notes.txt"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--outputs=raw.mrd"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--unknown"])
    end
end

@testitem "KomaUI" tags=[:koma] begin

    using Blink
    is_CI = Base.get_bool_env("CI", false)

    @static if Sys.islinux()
        KomaMRI.enable_unsafe_electron(!is_CI)
    end

    # Unfortunately Blink does not work on macOS in GitHub's CI
    # https://github.com/JuliaGizmos/Blink.jl/issues/325
    if !(Sys.isapple() && is_CI)
    # Opens UI
    w = KomaUI(return_window=true)
    @testset "Open UI" begin
        @test "index" == @js w document.getElementById("content").dataset.content
    end

    @testset "PulsesGUI" begin
        @js w document.getElementById("button_pulses_seq").click()
        @test "sequence" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_pulses_kspace").click()
        @test "kspace" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_pulses_M0").click()
        @test "m0" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_pulses_M1").click()
        @test "m1" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_pulses_M2").click()
        @test "m2" == @js w document.getElementById("content").dataset.content
    end

    @testset "PhantomGUI" begin
        @js w document.getElementById("button_phantom").click()
        @test "phantom" == @js w document.getElementById("content").dataset.content
    end

    @testset "ParamsGUI" begin
        @js w document.getElementById("button_scanner").click()
        @test "scanneparams" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_sim_params").click()
        @test "simparams" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_rec_params").click()
        @test "recparams" == @js w document.getElementById("content").dataset.content
    end

    @testset "Simulation" begin
        @js w document.getElementById("simulate!").click()
        @test "sig" == @js w document.getElementById("content").dataset.content
    end

    @testset "SignalGUI" begin
        @js w document.getElementById("button_sig").click()
        @test "sig" == @js w document.getElementById("content").dataset.content
    end

    @testset "Reconstruction" begin
        @js w document.getElementById("recon!").click()
        @test "absi" == @js w document.getElementById("content").dataset.content
    end

    @testset "ReconGUI" begin
        @js w document.getElementById("button_reconstruction_absI").click()
        @test "absi" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_reconstruction_angI").click()
        @test "angi" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_reconstruction_absK").click()
        @test "absk" == @js w document.getElementById("content").dataset.content
    end

    @testset "ExportToMAT" begin
        @js w document.getElementById("button_matfolder").click()
        @test "matfolder" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_matfolderseq").click()
        @test "matfolderseq" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_matfolderpha").click()
        @test "matfolderpha" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_matfoldersca").click()
        @test "matfoldersca" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_matfolderraw").click()
        @test "matfolderraw" == @js w document.getElementById("content").dataset.content

        @js w document.getElementById("button_matfolderima").click()
        @test "matfolderima" == @js w document.getElementById("content").dataset.content
    end

    if !isnothing(w)
        close(w)
    end

    end # if !(Sys.isapple() && is_CI)
end
