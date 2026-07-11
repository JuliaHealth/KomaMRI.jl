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

        opts = KomaMRI.parse_cli_args(["-i", "scanner.sys", "brain.h5", "epi.seq", "-o", "image.mat", "raw.mrd"])
        @test opts.sequence == "epi.seq"
        @test opts.phantom == "brain.h5"
        @test opts.scanner == "scanner.sys"
        @test opts.sim_output == "raw.mrd"
        @test opts.recon_output == "image.mat"

        opts = KomaMRI.parse_cli_args(["-i", "_", "epi.seq"])
        @test opts.sequence == "epi.seq"
        @test isnothing(opts.phantom)

        opts = KomaMRI.parse_cli_args(["-o", "raw.mat"])
        @test opts.sim_output == "raw.mat"
        @test isnothing(opts.recon_output)

        opts = KomaMRI.parse_cli_args(["-o", "_"])
        @test isnothing(opts.sim_output)
        @test isnothing(opts.recon_output)
    end

    @testset "Parameters and backend" begin
        opts = KomaMRI.parse_cli_args([
            "--backend=CPU",
            "-s", "gpu=true",
            "-s", "Nthreads=4",
            "-s", "sim_method=BlochMagnus4",
            "-s", "precision=f32",
            "-s", "max_rf_block_length=Inf",
            "-s", "offsets=[1,2.5,true,label]",
            "-r", "reco=direct",
            "-r", "shape=(2,3)",
        ])
        @test opts.backend == "CPU"
        @test opts.sim_params["gpu"] == true
        @test opts.sim_params["Nthreads"] == 4
        @test opts.sim_params["sim_method"] isa KomaMRI.BlochMagnus4
        @test opts.sim_params["precision"] == "f32"
        @test opts.sim_params["max_rf_block_length"] == Inf
        @test opts.sim_params["offsets"] == [1, 2.5, true, "label"]
        @test opts.recon_params[:reco] == "direct"
        @test opts.recon_params[:shape] == (2, 3)

        KomaMRI.load_cli_backend!(opts)
        @test opts.sim_params["gpu"] == false

        @test_throws ErrorException KomaMRI.load_cli_backend!(KomaMRI.CLIOptions(backend="NoBackend"))
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

    @testset "Input and output files" begin
        path = @__DIR__
        repo = dirname(path)
        phantom_file = joinpath(repo, "KomaMRIFiles", "test", "test_files", "phantom", "brain_nomotion_w.phantom")
        jemris_file = joinpath(repo, "KomaMRIFiles", "test", "test_files", "phantom", "column1d.h5")
        @test KomaMRI.load_cli_phantom(phantom_file) isa KomaMRI.Phantom
        @test KomaMRI.load_cli_phantom(jemris_file) isa KomaMRI.Phantom

        sys, seq, obj = KomaMRI.cli_inputs(KomaMRI.CLIOptions(scanner="scanner.sys"))
        @test sys isa KomaMRI.Scanner
        @test seq isa KomaMRI.Sequence
        @test obj isa KomaMRI.Phantom

        raw = RawAcquisitionData(ISMRMRDFile(joinpath(path, "test_files", "Koma_signal.mrd")))
        dir = mktempdir()
        @test KomaMRI.cli_output_dir("raw.mrd") == "."
        @test KomaMRI.mk_cli_output_dir(joinpath(dir, "nested", "raw.mrd")) == joinpath(dir, "nested")
        @test isdir(joinpath(dir, "nested"))

        raw_mrd = joinpath(dir, "raw.mrd")
        raw_mat = joinpath(dir, "raw.mat")
        img_mat = joinpath(dir, "image.mat")
        KomaMRI.save_cli_raw(raw, raw_mrd)
        KomaMRI.save_cli_raw(raw, raw_mat)
        KomaMRI.save_cli_recon(rand(2, 2, 1), Dict{Symbol,Any}(:reco => "direct"), img_mat)
        @test isfile(raw_mrd)
        @test isfile(raw_mat)
        @test isfile(img_mat)
    end

    @testset "Batch execution" begin
        dir = mktempdir()
        raw_mrd = joinpath(dir, "raw.mrd")
        img_mat = joinpath(dir, "image.mat")
        KomaMRI.run_cli(KomaMRI.CLIOptions(sim_output=raw_mrd, recon_output=img_mat))
        @test isfile(raw_mrd)
        @test isfile(img_mat)
    end

    @testset "App help" begin
        @static if VERSION >= v"1.12"
            @test occursin("KomaMRI command line app.", KomaMRI.CLI_HELP)
            redirect_stdout(devnull) do
                KomaMRI.CLI.print_help()
                @test isnothing(KomaMRI.CLI.main(["--help"]))
            end
            @test_throws ErrorException KomaMRI.CLI.main(["--unknown"])
        end

        @test isnothing(KomaMRI.print_cli_versions())
        @test KomaMRI.load_cli_preferences!(KomaMRI.CLIOptions()) isa KomaMRI.CLIOptions
    end

    @testset "Errors" begin
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs", "a.seq", "b.phantom", "c.sys", "d.seq"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--inputs", "notes.txt"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--outputs=raw.mrd"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["--unknown"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["-s", "sim_method=NoMethod"])
        @test_throws ErrorException KomaMRI.parse_cli_args(["-s", "sim_method=Phantom"])
        @test_throws ErrorException KomaMRI.run_cli(["--unknown"])
        @test_throws ErrorException KomaMRI.run_cli(["-b", "NoBackend"])
        @test_throws ErrorException KomaMRI.load_cli_phantom("brain.txt")
        @test_throws ErrorException KomaMRI.save_cli_raw(nothing, "raw.txt")
        @test_throws ErrorException KomaMRI.save_cli_recon(nothing, Dict{Symbol,Any}(), "image.nii")
    end
end

@testitem "KomaUI" tags=[:koma] begin
    @testset "MAT export" begin
        scanner = Scanner()
        sequence = KomaMRI.setup_sequence(scanner)
        phantom = KomaMRI.setup_phantom()
        raw = KomaMRI.setup_raw()
        image = ComplexF64[0 0; 0 0]
        rec_params = Dict{Symbol,Any}(:reco => "direct")

        dir = mktempdir()
        message = KomaMRI.export_2_mat(sequence, phantom, scanner, raw, rec_params, image, dir; type="image")
        @test readdir(dir) == ["data_image.mat"]
        @test occursin("<b>Name:</b> data_image.mat", message)

        dir = mktempdir()
        message = KomaMRI.export_2_mat(sequence, phantom, scanner, raw, rec_params, image, dir; type="sequence")
        @test sort(readdir(dir)) == ["data_kspace.mat", "data_moments.mat", "data_sequence.mat"]
        @test occursin("<b>Names:</b> data_sequence.mat, data_kspace.mat, data_moments.mat", message)
    end

    @testset "Rendered desktop UI" begin
        is_CI = Base.get_bool_env("CI", false)
        if Sys.isapple() && is_CI
            @test_skip "Electron windows are unavailable on macOS CI"
        else
            w = KomaUI(; return_window=true, sim=Dict{String,Any}("gpu" => false))
            ui = KomaMRI
            session = w.session[]
            code(source) = ui.Bonito.JSCode(source)
            value(source) = ui.Bonito.evaljs_value(session, code(source))
            function wait_until(test; timeout=10)
                start = time()
                while !test()
                    time() - start > timeout && return false
                    sleep(0.05)
                end
                return true
            end

            try
                @test value("document.querySelectorAll('.koma-nav-icon').length") >= 6
                @test value("getComputedStyle(document.querySelector('.koma-nav-icon'), '::before').content !== 'none'")
                @test value("document.querySelector('#seqfilepicker input[type=file]') !== null")
                @test value("document.querySelector('#seqfilepicker .koma-file-name').textContent") == "epi.seq"

                upload_name = "label_test_with_a_deliberately_long_filename.seq"
                upload_data = read(
                    joinpath(pkgdir(KomaMRI), "KomaMRIFiles", "test", "test_files", "pulseq", "basic_tests", "v1.4", "label_test.seq"),
                )
                data = join(upload_data, ',')
                ui.Bonito.evaljs(session, code("""
                    (() => {
                        const input = document.querySelector('#seqfilepicker input[type=file]');
                        const file = new File([new Uint8Array([$data])], $(repr(upload_name)));
                        Object.defineProperty(input, 'files', {value: [file], configurable: true});
                        input.dispatchEvent(new Event('change', {bubbles: true}));
                    })()
                """))
                @test wait_until(() -> value("document.querySelector('#seqfilepicker .koma-file-name').title") == upload_name)
                @test value("document.querySelector('#seqfilepicker .koma-file-name').textContent") == KomaMRI.display_filename(upload_name)
                @test wait_until(() -> value("document.getElementById('toastTitle0').textContent.includes('...')"))
                @test wait_until(() -> w.state[] == "sequence")
                @test wait_until(() -> value("document.querySelector('.js-plotly-plot') !== null"))

                value("document.getElementById('button_pulses_seq').click(); true")
                @test wait_until(() -> w.state[] == "sequence")
                @test wait_until(() -> value("document.querySelector('.js-plotly-plot') !== null"))
                @test value("(() => { const plot = document.querySelector('.js-plotly-plot'); return Math.abs(plot.clientWidth - plot.parentElement.clientWidth) <= 1; })()")

                ui.toast!(w, 0, "Loaded <b>epi.seq</b>", "Ready to <b>simulate</b>?")
                @test wait_until(() -> value("document.getElementById('toastTitle0').innerHTML.includes('epi.seq')"))
                @test value("document.getElementById('toastBody0').innerHTML.includes('simulate')")

                ui.start_simulation_progress!(w)
                @test value("document.getElementById('simulate!').disabled")
                @test value("document.getElementById('simul_progress') !== null")
                ui.finish_simulation_progress!(w)
                @test !value("document.getElementById('simulate!').disabled")

                seq_ui[] = PulseDesigner.EPI_example(; sys=sys_ui[])
                @test wait_until(() -> w.state[] == "sequence")
                @test wait_until(() -> value("document.querySelector('.js-plotly-plot') !== null"))

                img_ui[] = reshape(ComplexF32[1, 0, 0, 1], 2, 2, 1)
                @test wait_until(() -> w.state[] == "absi")
                @test wait_until(() -> value("document.querySelector('.js-plotly-plot') !== null"))
            finally
                close(w)
            end
        end
    end
end
