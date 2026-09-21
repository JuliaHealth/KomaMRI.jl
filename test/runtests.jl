using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:koma in ti.tags) #verbose=true

# include("../KomaMRICore/test/runtests.jl")
# include("../KomaMRIPlots/test/runtests.jl")

@testitem "MRIReco recon" tags=[:koma] begin
    using MAT

    mat_values(value::Number) = [value]
    mat_values(value::AbstractString) = [value]
    mat_values(value) = vec(collect(value))

    function test_reconstruction_mat(result)
        mktempdir() do dir
            KomaMRI.export_2_mat_image(result, Dict(:reco=>"direct"), dir)
            saved = matread(joinpath(dir, "image.mat"))["reconstruction"]
            # MAT structs preserve each labeled complex array, its dimensions, and reconstruction policy.
            @test saved["policy"] == Dict(string(k)=>string(v) for (k, v) in pairs(result.policy)) &&
                length(saved["images"]) == length(result.images)
            for (stored, entry) in zip(saved["images"], result.images)
                @test mat_values(stored["size"]) == collect(size(entry.image)) &&
                    ntuple(d -> size(stored["data"], d), 6) == size(entry.image) &&
                    vec(stored["data"]) == vec(parent(entry.image))
                @test stored["labels"] == Dict(string(k)=>(k == :ROLE ? string(v) : v) for (k, v) in pairs(entry.labels)) &&
                    all(mat_values(stored["source"][string(k)]) == collect(v) for (k, v) in pairs(entry.source))
            end
        end
    end

    function test_raw_mat(raw, roles)
        mktempdir() do dir
            KomaMRI.export_2_mat_raw(raw, dir)
            saved = matread(joinpath(dir, "raw.mat"))["raw"]
            # Every coil sample, trajectory coordinate, and MRD header field survives the MAT roundtrip.
            @test mat_values(saved["params"]["encodedSize"]) == raw.params["encodedSize"] &&
                haskey(saved["params"], "userParameters") == haskey(raw.params, "userParameters")
            @test length(saved["profiles"]) == length(raw.profiles)
            @test all(stored["data"] == profile.data && stored["traj"] == profile.traj && stored["role"] == role &&
                all(mat_values(stored["head"][string(k)]) == mat_values(getproperty(profile.head, k))
                    for k in propertynames(profile.head) if k != :idx) &&
                all(mat_values(stored["head"]["idx"][string(k)]) == mat_values(getproperty(profile.head.idx, k))
                    for k in propertynames(profile.head.idx))
                for (stored, profile, role) in zip(saved["profiles"], raw.profiles, roles))
        end
    end

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

    @testset "Cartesian point-source export and reconstruction" begin
        fov = 0.08
        function cartesian_sequence(matrix; fov, duration=1e-3)
            nx, ny, nz = matrix
            scale = 1 / (γ * fov * duration)
            seq = Sequence()
            @addblock for partition in 0:(nz - 1), line in 0:(ny - 1)
                start = [-fld(nx, 2), line - fld(ny, 2), partition - fld(nz, 2)]
                pre = [Grad(k * scale, duration) for k in start]
                seq += (x=pre[1], y=pre[2], z=pre[3])
                seq += (x=Grad((nx - 1) * scale, duration), ADC(nx, duration),
                    LabelSet(line, "LIN"), nz > 1 ? LabelSet(partition, "PAR") : nothing)
                seq += (x=Grad(-(start[1] + nx - 1) * scale, duration), y=-pre[2], z=-pre[3])
            end
            seq.DEF = Dict("Nx"=>nx, "Ny"=>ny, "Nz"=>nz, "FOV"=>fill(fov, 3))
            return seq
        end
        function reconstruct_point(seq, position)
            _, k = get_kspace(seq)
            signal = reshape(ComplexF32.(exp.(-2π * im .* (k * position))), :, 1)
            raw = signal_to_raw_data(signal, seq)
            return reconstruct_with_labels(raw; rec_params=Dict(:reco=>"direct"))
        end

        # A point three voxels from isocenter remains there for x/y/z readouts, with no padded rows.
        samples, displacement = 16, 3
        rotations = ([1 0 0; 0 1 0; 0 0 1], [0 0 1; 1 0 0; 0 1 0], [0 1 0; 0 0 1; 1 0 0])
        for (axis, rotation) in enumerate(rotations)
            seq = rotation * cartesian_sequence((samples, 1, 1); fov)
            for (dimension, name) in enumerate(("Nx", "Ny", "Nz"))
                seq.DEF[name] = dimension == axis ? samples : 1
            end
            image = only(reconstruct_point(seq, rotation * [displacement * fov / samples, 0, 0]).images).image
            @test size(image)[1:3] == (samples, 1, 1)
            @test argmax(abs.(image)) == CartesianIndex(samples ÷ 2 + 1 + displacement, 1, 1, 1, 1, 1)
        end

        # Physical 3D encoding must locate an off-center point without PAR labels.
        matrix, displacement = (8, 8, 8), [1, -2, 2]
        seq = cartesian_sequence(matrix; fov)
        position = displacement .* (fov ./ collect(matrix))
        for extensions in seq.EXT
            filter!(ext -> !(ext isa LabelSet && ext.labelstring == "PAR"), extensions)
        end
        volume = reconstruct_point(seq, position)
        unlabeled = only(volume.images).image
        @test argmax(abs.(unlabeled)) == CartesianIndex((matrix .÷ 2 .+ 1 .+ Tuple(displacement))..., 1, 1, 1)
        test_reconstruction_mat(volume)

        # A varying spatial direction without a matrix must fail rather than invent an image size.
        missing_matrix = signal_to_raw_data(ones(ComplexF32, prod(matrix), 1), seq)
        missing_matrix.params["encodedSize"][3] = 1
        @test_throws ArgumentError reconstruct_with_labels(missing_matrix)

        @testset "Mixed 2D slices and 1D navigators" begin
            image_samples = 8
            slice_seq = cartesian_sequence((image_samples, image_samples, 1); fov)
            imaging_seq = Sequence()
            @addblock for slice in 0:1
                imaging_seq += (LabelSet(slice, "SLC"),)
                imaging_seq += slice_seq
            end
            imaging_seq.DEF = copy(slice_seq.DEF)
            point_position = (fov / image_samples) .* [1, -1, 1]
            _, imaging_k = get_kspace(imaging_seq)
            imaging_signal = reshape(ComplexF32.(exp.(-2π * im .* (imaging_k * point_position))), :, 1)
            imaging_signal[(image_samples^2 + 1):end, :] .*= 2
            imaging_raw = signal_to_raw_data(imaging_signal, imaging_seq)
            rec_params = Dict(:reco=>"direct")
            reference = reconstruct_with_labels(imaging_raw; rec_params).images

            # A wider z navigator must not rescale either slice or turn the images into a 3D volume.
            # Export normalization differs by a factor of three, including Float32 coordinate rounding.
            navigator = last(rotations) * cartesian_sequence((image_samples, 1, 1); fov=fov / 3)
            mixed_seq = copy(imaging_seq)
            @addblock begin
                mixed_seq += (LabelSet(1, "NAV"), LabelSet(0, "SLC"))
                mixed_seq += navigator
            end
            mixed_seq.DEF = copy(imaging_seq.DEF)
            _, mixed_k = get_kspace(mixed_seq)
            navigator_k = @view mixed_k[(size(imaging_k, 1) + 1):end, :]
            navigator_signal = reshape(ComplexF32.(exp.(-2π * im .* (navigator_k * point_position))), :, 1)
            mixed_raw = signal_to_raw_data(vcat(imaging_signal, navigator_signal), mixed_seq)
            result = reconstruct_with_labels(mixed_raw; rec_params)
            images = filter(entry -> entry.labels.ROLE === :imaging, result.images)
            line = only(entry for entry in result.images if entry.labels.ROLE === :navigator).image
            @test [entry.labels.SLC for entry in images] == [0, 1]
            @test all(parent(image.image) ≈ parent(original.image) for (image, original) in zip(images, reference))
            @test parent(images[2].image) ≈ 2 .* parent(images[1].image)
            @test size(line)[1:3] == (image_samples, 1, 1)
            navigator_displacement = 3
            @test argmax(abs.(line)) == CartesianIndex(image_samples ÷ 2 + 1 + navigator_displacement, 1, 1, 1, 1, 1)

            test_reconstruction_mat(result)
            foreign_raw = RawAcquisitionData(copy(mixed_raw.params), mixed_raw.profiles)
            delete!(foreign_raw.params, "userParameters")
            test_raw_mat(foreign_raw, [fill("imaging", length(imaging_raw.profiles)); "navigator"])
            mktempdir() do dir
                KomaMRI.export_2_mat_sequence(mixed_seq, dir)
                saved = matread(joinpath(dir, "seq_sequence.mat"))["sequence"]
                context = adc_label_context(mixed_seq)
                # Sequence export retains accumulated NAV/SLC/LIN labels aligned with ADC blocks and samples.
                @test mat_values(saved["adc"]["block"]) == context.adc_blocks &&
                    mat_values(saved["adc"]["num_samples"]) == mixed_seq.ADC[context.adc_blocks].N &&
                    all(mat_values(saved["adc"]["labels"][string(k)]) == getproperty.(context.labels[context.adc_blocks], k)
                        for k in propertynames(first(context.labels)))
                @test all(mat_values(saved["definitions"][k]) == mat_values(v) for (k, v) in mixed_seq.DEF)
            end

            # Discarded off-axis samples outside normalized ±0.5 cannot change geometry or scaling.
            padded_raw = deepcopy(mixed_raw)
            padded_nav = last(padded_raw.profiles)
            discarded_coordinate = Float32[0, 1, 1]
            discarded_signal = zeros(ComplexF32, 1, size(padded_nav.data, 2))
            padded_nav.traj = hcat(-discarded_coordinate, padded_nav.traj, discarded_coordinate)
            padded_nav.data = vcat(discarded_signal, padded_nav.data, discarded_signal)
            padded_nav.head.discard_pre = padded_nav.head.discard_post = 1
            padded_nav.head.number_of_samples += 2
            padded_nav.head.center_sample += 1
            padded = reconstruct_with_labels(padded_raw; rec_params).images
            @test parent(only(entry for entry in padded if entry.labels.ROLE === :navigator).image) ≈ parent(line)
        end
    end

    @testset "Shared reconstruction display ranges" begin
        data = ones(ComplexF32, 2, 2, 2, 2, 2, 1)
        data[:, :, 2, :, :, :] .*= 2
        data[:, :, :, 2, :, :] .*= 10
        data[:, :, :, :, 2, :] .*= 3
        images = [
            (; labels=(; ROLE=role, DIM=2, SET=set, AVG=average),
               source=(; ECO=[0, 5]), image=amplitude .* data)
            for (role, set, average, amplitude) in (
                (:imaging, 0, 0, 1), (:imaging, 0, 1, 4),
                (:imaging, 1, 0, 100), (:navigator, 0, 0, 1000),
            )
        ]
        cached = KomaMRI._with_magnitude_limits(images)

        # Average, partition, and coil changes share a range; ECO, SET, and ROLE do not.
        @test cached[1].magnitude_limits == cached[2].magnitude_limits ==
            [(4, 96), (40, 960)]
        @test cached[3].magnitude_limits == [(400, 2400), (4000, 24000)]
        @test cached[4].magnitude_limits == [(4000, 24000), (40000, 240000)]
    end

    @testset "Reconstruction policy" begin
        average_profiles = map(raw.profiles) do profile
            head = deepcopy(profile.head)
            head.idx.average = 1
            Profile(head, profile.traj, 3im .* profile.data)
        end
        average_raw = RawAcquisitionData(raw.params, [raw.profiles; average_profiles])
        rec_params = Dict{Symbol,Any}(:reco => "direct")
        separate = KomaMRI.reconstruct_with_labels(average_raw; rec_params)
        combined = KomaMRI.reconstruct_with_labels(
            average_raw; recon_policy=(; AVG=:combine), rec_params
        )
        first_average = only(image for image in separate.images if image.labels.AVG == 0)
        second_average = only(image for image in separate.images if image.labels.AVG == 1)
        reference = parent(first_average.image)

        # AVG separation preserves phase; combining sums complex signals without normalization.
        @test parent(second_average.image) ≈ 3im .* reference
        @test parent(only(combined.images).image) ≈ (1 + 3im) .* reference

        coil_profiles = map(raw.profiles) do profile
            head = deepcopy(profile.head)
            head.available_channels = 2
            head.active_channels = 2
            Profile(head, profile.traj, hcat(3im .* profile.data, 4im .* profile.data))
        end
        coil_raw = RawAcquisitionData(raw.params, coil_profiles)
        coil_result = KomaMRI.reconstruct_with_labels(coil_raw; rec_params)
        rss_result = KomaMRI.reconstruct_with_labels(
            coil_raw; recon_policy=(; COIL=:rss), rec_params
        )
        coil_axis, coil_rss = only(coil_result.images).image, only(rss_result.images).image

        # Coil gains 3im and 4im retain phase separately: RSS gives magnitude 5, coherent sum 7.
        @test parent(coil_axis)[:, :, :, :, 1:1, :] ≈ 3im .* reference
        @test parent(coil_axis)[:, :, :, :, 2:2, :] ≈ 4im .* reference
        @test parent(coil_rss) ≈ 5 .* abs.(reference)
        test_reconstruction_mat(coil_result)
        test_reconstruction_mat(rss_result)
        test_raw_mat(coil_raw, fill("imaging", length(coil_raw.profiles)))
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

        # A scanner input replaces the defaults, including receive channels and hardware limits.
        scanner_file = joinpath(mktempdir(), "scanner.sys")
        expected = Scanner(limits=HardwareLimits(B0=3.0), receiver=BirdcageCoilSens(ncoils=4))
        write_scanner(expected, scanner_file)
        sys, seq, obj = KomaMRI.cli_inputs(KomaMRI.CLIOptions(scanner=scanner_file))
        @test sys.limits.B0 == expected.limits.B0
        @test get_sens(sys.receiver, [0.0], [0.0], [0.0]) ==
            get_sens(expected.receiver, [0.0], [0.0], [0.0])

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
        scanner_file = joinpath(dir, "scanner.sys")
        maps = cat(fill(1.0 + 0im, 2, 2, 2), fill(0.0 + 2im, 2, 2, 2); dims=4)
        write_scanner(Scanner(receiver=ArbitraryCoilSens([-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], maps)), scanner_file)
        KomaMRI.run_cli(KomaMRI.CLIOptions(scanner=scanner_file, sim_output=raw_mrd, recon_output=img_mat))
        # The CLI simulation uses the saved coil maps, including their relative complex gain.
        raw = RawAcquisitionData(ISMRMRDFile(raw_mrd))
        @test all(profile.data[:, 2] == 2im .* profile.data[:, 1] for profile in raw.profiles)
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

@testitem "Scanner MAT export" tags=[:koma] begin
    using MAT

    coordinates = range(-0.1, 0.1; length=2)
    maps = reshape(ComplexF32.(1:16) .* (1 + 2im), 2, 2, 2, 2)
    sys = Scanner(limits=HardwareLimits(B0=3.0),
        receiver=ArbitraryCoilSens(coordinates, coordinates, coordinates, maps))
    mktempdir() do dir
        KomaMRI.export_2_mat_scanner(sys, dir)
        saved = matread(joinpath(dir, "scanner.mat"))["scanner"]
        # MATLAB receives all four scanner components, with SI limits and complex maps intact.
        @test Set(keys(saved)) == Set(["limits", "gradient", "receiver", "transmitter"])
        @test all(saved["limits"][replace(string(name), "Δ" => "d")] == getproperty(sys.limits, name)
            for name in fieldnames(HardwareLimits))
        @test saved["gradient"]["type"] == "LinearXYZ" && saved["transmitter"]["type"] == "UniformTransmit"
        @test saved["receiver"]["type"] == "ArbitraryCoilSens"
        @test saved["receiver"]["coil_sens"] == maps
        @test all(vec(saved["receiver"][axis]) == coordinates for axis in ("x", "y", "z"))
    end
end

@testitem "KomaUI" tags=[:koma] begin
    using Bonito

    triggered = Sequence()
    @addblock triggered += PulseDesigner.make_trigger(:physio1; duration=1e-3)

    @testset "MAT export" begin
        scanner = Scanner()
        sequence = KomaMRI.setup_sequence(scanner)
        phantom = KomaMRI.setup_phantom()
        raw = KomaMRI.setup_raw()
        image = ComplexF64[1 im; -2 3im]
        rec_params = Dict{Symbol,Any}(:reco => "direct")

        dir = mktempdir()
        message = KomaMRI.export_2_mat(sequence, phantom, scanner, raw, rec_params, image, dir; type="image")
        @test readdir(dir) == ["data_image.mat"]
        @test occursin("<b>Name:</b> data_image.mat", message)
        # Plain-array callers retain the legacy image variable and complex data.
        @test KomaMRI.matread(joinpath(dir, "data_image.mat"))["image"] == image

        dir = mktempdir()
        message = KomaMRI.export_2_mat(sequence, phantom, scanner, raw, rec_params, image, dir; type="sequence")
        @test sort(readdir(dir)) == ["data_kspace.mat", "data_moments.mat", "data_sequence.mat"]
        @test occursin("<b>Names:</b> data_sequence.mat, data_kspace.mat, data_moments.mat", message)
    end

    @testset "Window-owned data" begin
        w = KomaUI(; seq=triggered, obj=Phantom(x=[0.0]), show_window=false,
            return_window=true, verbose=false, sim=Dict("gpu" => false))
        other = KomaUI(; seq=Sequence(), obj=Phantom(x=[0.0]), show_window=false,
            return_window=true, verbose=false, sim=Dict("gpu" => false))
        try
            # Opening another window must not reset the first window's acquisition or physiology.
            @test w.seq[] === triggered && w.physio[].period == 1.0

            # Every data observable and its view update are isolated to their owning window.
            updates = (
                seq=triggered,
                obj=Phantom(x=[0.0], ρ=[2.0]),
                sys=Scanner(limits=HardwareLimits(B0=3.0)),
                physio=CardiacSignal(; heart_rate=1.25),
                raw=KomaMRI.setup_raw(),
                img=ComplexF32[0 1; 2 3],
                sim_params=Dict("gpu" => false, "precision" => "f64"),
                rec_params=Dict(:reco => "standard", :iterations => 0),
            )
            for (field, value) in pairs(updates)
                other_value = getproperty(other, field)[]
                other_content = other.content[]
                getproperty(w, field)[] = value
                @test getproperty(other, field)[] === other_value && other.content[] === other_content
            end

            # File loading changes neither the other window's acquisition nor its reload path.
            other_sequence = other.seq[]
            load_file!(w, :sequence, joinpath(pkgdir(KomaMRI), "KomaMRIFiles", "test",
                "test_files", "pulseq", "basic_tests", "v1.4", "label_test.seq"))
            @test other.seq[] === other_sequence && isempty(other.files)

            # Closing one window must leave the other's observable subscriptions active.
            close(w)
            other.seq[] = triggered
            @test other.state[] == "sequence" && other.physio[].period == 1.0
        finally
            close(w)
            close(other)
        end
    end

    @testset "Rendered desktop UI" begin
        is_CI = Base.get_bool_env("CI", false)
        if Sys.isapple() && is_CI
            @test_skip "Electron windows are unavailable on macOS CI"
        else
            w = KomaUI(; return_window=true, sim=Dict{String,Any}("gpu" => false))
            session = w.session[]
            click_button(id) =
                Bonito.evaljs(session, js"document.getElementById($(id)).click()")
            plot_rendered(state) = Bonito.evaljs_value(
                session,
                js"""
                    document.getElementById('content').dataset.content === $(state) &&
                        document.querySelector('#content .js-plotly-plot') !== null
                """,
            )
            range_slider_visible() = Bonito.evaljs_value(
                session,
                js"""
                    document.querySelector('#content .rangeslider-container')
                        .getBoundingClientRect().bottom <=
                    document.getElementById('content').getBoundingClientRect().bottom
                """,
            )
            try
                @testset "Open UI" begin
                    @test w.state[] == "index"
                    # Both expanded and collapsed sidebar logos must resolve to loaded images.
                    @test timedwait(() -> Bonito.evaljs_value(session, js"""
                        ['.koma-logo-full', '.koma-logo-spin'].every(selector => {
                            const image = document.querySelector(selector);
                            return image !== null && image.complete && image.naturalWidth > 0;
                        })
                    """), 30) == :ok
                end

                @testset "Julia actions" begin
                    # Julia actions select the same views as the menu, without a DOM click.
                    click!(w, :view_hardware_limits)
                    @test w.state[] == "scanneparams"
                    click!(w, :home)
                    @test w.state[] == "index"
                    @test_throws KeyError click!(w, :unknown_action)
                end

                @testset "Sequence views" begin
                    click_button("button_pulses_seq")
                    @test timedwait(() -> w.state[] == "sequence", 30) == :ok
                    @test timedwait(() -> plot_rendered("sequence"), 30) == :ok
                    @test range_slider_visible()
                    @test Bonito.evaljs_value(
                        session, js"document.getElementById('main').clientHeight === window.innerHeight"
                    )

                    click_button("button_pulses_kspace")
                    @test timedwait(() -> w.state[] == "kspace", 30) == :ok
                    @test timedwait(() -> plot_rendered("kspace"), 30) == :ok

                    # Moments and slew rate share one time slider with the sequence panel above.
                    for (button, state) in (("M0", "m0"), ("M1", "m1"),
                        ("M2", "m2"), ("slew_rate", "slew_rate"))
                        click_button("button_pulses_$button")
                        @test timedwait(() -> w.state[] == state, 30) == :ok
                        @test timedwait(() -> Bonito.evaljs_value(session, js"""
                            (() => {
                                const panels = document.querySelectorAll('#content .koma-time-plot');
                                return panels.length === 2 &&
                                    panels[0].querySelector('.legend') !== null &&
                                    document.querySelectorAll('#content .rangeslider-container').length === 1;
                            })()
                        """), 30) == :ok
                        @test range_slider_visible()
                    end
                end

                @testset "Phantom and parameters" begin
                    click_button("button_phantom")
                    @test timedwait(() -> w.state[] == "phantom", 30) == :ok

                    click_button("button_scanner")
                    @test timedwait(() -> w.state[] == "coils", 30) == :ok
                    @test timedwait(() -> plot_rendered("coils"), 30) == :ok

                    click_button("button_scanner_params")
                    @test timedwait(() -> w.state[] == "scanneparams", 30) == :ok

                    click_button("button_sim_params")
                    @test timedwait(() -> w.state[] == "simparams", 30) == :ok

                    click_button("button_rec_params")
                    @test timedwait(() -> w.state[] == "recparams", 30) == :ok
                end

                @testset "Simulation and raw signal" begin
                    # Simulation must use replacement options, not the launch-time dictionary.
                    w.sim_params[] = Dict("gpu" => false, "precision" => "f64")
                    @test w.state[] == "simparams"
                    click_button("simulate!")
                    @test timedwait(() -> w.state[] == "sig", 180) == :ok
                    @test timedwait(() -> plot_rendered("sig"), 30) == :ok
                    @test !isempty(w.raw[].profiles)
                    @test w.raw[].params["userParameters"]["precision"] == "f64"

                    click_button("button_scanner_params")
                    @test timedwait(() -> w.state[] == "scanneparams", 30) == :ok

                    click_button("button_sig")
                    @test timedwait(() -> w.state[] == "sig", 30) == :ok
                    @test timedwait(() -> plot_rendered("sig"), 30) == :ok
                    @test range_slider_visible()

                    # Reload restores the saved simulation samples and reopens the raw-data plot.
                    simulated_samples = copy(first(w.raw[].profiles).data)
                    fill!(first(w.raw[].profiles).data, 0)
                    click_button("button_scanner_params")
                    @test timedwait(() -> w.state[] == "scanneparams", 30) == :ok
                    click_button("button_reload_raw")
                    @test timedwait(() -> first(w.raw[].profiles).data == simulated_samples, 30) == :ok
                    @test timedwait(() -> plot_rendered("sig"), 30) == :ok
                end

                @testset "Reconstruction and image views" begin
                    # Zero iterations retains the zero initial image, proving replacement options are used.
                    w.rec_params[] = Dict(:reco => "standard", :iterations => 0)
                    @test w.state[] == "recparams"
                    click_button("recon!")
                    @test timedwait(() -> w.state[] == "absi", 180) == :ok
                    @test all(iszero, first(w.img[].images).image)

                    # Switching back to direct reconstruction restores an image from the nonzero signal.
                    w.rec_params[] = Dict(:reco => "direct")
                    click_button("recon!")
                    @test timedwait(() -> w.state[] == "absi", 180) == :ok
                    @test timedwait(() -> plot_rendered("absi"), 30) == :ok
                    @test !isempty(w.img[])
                    @test any(!iszero, first(w.img[].images).image)

                    # Components switch within the existing chart, without a separate UI page.
                    for (component, index) in (("Phase", 1), ("Magnitude", 0))
                        Bonito.evaljs(session, js"""
                            Array.from(document.querySelectorAll('#content .updatemenu-button'))
                                .find(button => button.textContent === $(component))
                                .dispatchEvent(new MouseEvent('click', {bubbles: true}));
                        """)
                        @test timedwait(30) do
                            Bonito.evaljs_value(session, js"""
                                document.querySelector('#content .js-plotly-plot')
                                    ._fullData[$(index)].visible === true
                            """)
                        end == :ok
                        @test w.state[] == "absi"
                    end

                    click_button("button_reconstruction_absK")
                    @test timedwait(() -> w.state[] == "absk", 30) == :ok
                    @test timedwait(() -> plot_rendered("absk"), 30) == :ok
                end

                @testset "Image export" begin
                    display = w.display[]
                    output = joinpath(tempdir(), "data_image.mat")
                    rm(output; force=true)
                    w.display[] = nothing
                    try
                        click_button("button_matfolderima")
                        @test timedwait(() -> w.state[] == "matfolderima", 30) == :ok
                        @test isfile(output)
                    finally
                        rm(output; force=true)
                        w.display[] = display
                    end
                end

                @testset "Observable updates" begin
                    w.seq[] = PulseDesigner.EPI_example(; sys=w.sys[])
                    @test timedwait(() -> w.state[] == "sequence", 30) == :ok
                    @test timedwait(() -> plot_rendered("sequence"), 30) == :ok

                    w.seq[] = triggered
                    @test timedwait(() -> w.state[] == "sequence", 30) == :ok
                    @test timedwait(() -> plot_rendered("sequence"), 30) == :ok
                    @test w.physio[].period == 1.0
                    @test Bonito.evaljs_value(
                        session,
                        js"document.querySelector('#content .rangeslider-container') === null",
                    )

                    w.physio[] = CardiacSignal(; heart_rate=1.25)
                    @test w.physio[].period == 0.8

                    w.seq[] = PulseDesigner.EPI_example(; sys=w.sys[])
                    @test w.physio[] == NoPhysioSignal()

                    w.obj[] = KomaMRI.setup_phantom()
                    @test timedwait(() -> w.state[] == "phantom", 30) == :ok
                    @test timedwait(() -> plot_rendered("phantom"), 30) == :ok

                    w.sys[] = Scanner()
                    @test timedwait(() -> w.state[] == "coils", 30) == :ok

                    # Scanner edits open the changed component, including notified nested limits.
                    w.sys[].limits.B0 = 3.0
                    notify(w.sys)
                    @test timedwait(() -> w.state[] == "scanneparams", 30) == :ok
                    sys = w.sys[]
                    w.sys[] = Scanner(; limits=sys.limits, gradient=sys.gradient,
                        transmitter=sys.transmitter, receiver=BirdcageCoilSens())
                    @test timedwait(() -> plot_rendered("coils"), 30) == :ok

                    w.raw[] = RawAcquisitionData(
                        ISMRMRDFile(joinpath(@__DIR__, "test_files", "Koma_signal.mrd"))
                    )
                    @test timedwait(() -> w.state[] == "sig", 30) == :ok
                    @test timedwait(() -> plot_rendered("sig"), 30) == :ok

                    # Slice changes reuse the graph; leaving the page releases its retained plot state.
                    w.img[] = cat(zeros(ComplexF32, 2, 2), ones(ComplexF32, 2, 2); dims=3)
                    @test timedwait(() -> w.state[] == "absi", 30) == :ok
                    @test timedwait(() -> plot_rendered("absi"), 30) == :ok
                    try
                        Bonito.evaljs(session, js"""
                            window.komaTestPlot = document.querySelector('#content .js-plotly-plot');
                            const slider = document.querySelector('#content input[type=range]');
                            slider.value = '2';
                            slider.dispatchEvent(new Event('input', {bubbles: true}));
                        """)
                        @test timedwait(() -> Bonito.evaljs_value(session, js"""
                            window.komaTestPlot === document.querySelector('#content .js-plotly-plot') &&
                            window.komaTestPlot.data[0].z.flat().every(value => value === 4)
                        """), 30) == :ok
                        click_button("button_scanner")
                        @test timedwait(() -> Bonito.evaljs_value(session, js"""
                            !window.komaTestPlot.isConnected &&
                            window.komaTestPlot.data === undefined &&
                            window.komaTestPlot._fullLayout === undefined &&
                            window.komaTestPlot._responsiveChartHandler === undefined
                        """), 30) == :ok
                    finally
                        Bonito.evaljs(session, js"delete window.komaTestPlot")
                    end
                end

                @testset "File loading" begin
                    files = joinpath(pkgdir(KomaMRI), "KomaMRIFiles", "test", "test_files")
                    sequence_file = joinpath(
                        files, "pulseq", "basic_tests", "v1.4", "label_test.seq"
                    )
                    load_file!(w, :sequence, sequence_file)
                    @test any(ext -> ext isa LabelInc, Iterators.flatten(w.seq[].EXT))
                    @test timedwait(() -> w.state[] == "sequence", 30) == :ok
                    # Loaded filenames stay discoverable when the sidebar labels are collapsed.
                    @test Bonito.evaljs_value(session,
                        js"document.getElementById('seq-dropdown').title === 'Sequence: label_test.seq'")

                    # A failed load must preserve both the input and its reload target.
                    loaded_sequence = w.seq[]
                    @test_throws SystemError load_file!(w, :sequence, joinpath(files, "missing.seq"))
                    @test w.seq[] === loaded_sequence && w.files[:sequence] == abspath(sequence_file)

                    # Reload must reread the selected path and display the changed sequence.
                    reload_source = joinpath(mktempdir(), "reload.seq")
                    cp(sequence_file, reload_source)
                    _, selected_file = KomaMRI.filepicker_selection(Dict(
                        "name" => "reload.seq",
                        "path" => reload_source,
                        "data" => read(sequence_file),
                    ))
                    load_file!(w, :sequence, selected_file)
                    click_button("button_phantom")
                    @test timedwait(() -> w.state[] == "phantom", 30) == :ok
                    cp(
                        joinpath(files, "pulseq", "basic_tests", "v1.4", "epi.seq"),
                        reload_source;
                        force=true,
                    )
                    previous_sequence = w.seq[]
                    click_button("button_reload_seq")
                    @test timedwait(() -> w.seq[] !== previous_sequence, 30) == :ok
                    @test !any(ext -> ext isa LabelInc, Iterators.flatten(w.seq[].EXT))
                    @test timedwait(() -> w.state[] == "sequence", 30) == :ok
                    @test timedwait(() -> plot_rendered("sequence"), 30) == :ok

                    @test Bonito.evaljs_value(session,
                        js"document.getElementById('seq-dropdown').title === 'Sequence: reload.seq'")
                    phantom_file = joinpath(files, "phantom", "column1d.h5")
                    previous_phantom = w.obj[]
                    load_file!(w, :phantom, phantom_file)
                    @test w.obj[].name == "column1d.h5"
                    @test timedwait(() -> w.state[] == "phantom", 30) == :ok

                    w.obj[] = previous_phantom
                    click_button("button_reload_phantom")
                    @test timedwait(() -> w.obj[].name == "column1d.h5", 30) == :ok
                    @test timedwait(() -> w.state[] == "phantom", 30) == :ok
                    @test timedwait(() -> plot_rendered("phantom"), 30) == :ok

                    # Scanner loading updates the receive view; reloading rereads changed limits and channels.
                    scanner_file = joinpath(mktempdir(), "scanner.sys")
                    write_scanner(Scanner(receiver=BirdcageCoilSens(ncoils=4)), scanner_file)
                    load_file!(w, :scanner, scanner_file)
                    @test get_n_coils(w.sys[].receiver) == 4
                    @test timedwait(() -> plot_rendered("coils"), 30) == :ok
                    write_scanner(Scanner(limits=HardwareLimits(B0=3.0)), scanner_file)
                    click_button("button_reload_scanner")
                    @test timedwait(() -> w.sys[].limits.B0 == 3.0 && get_n_coils(w.sys[].receiver) == 1, 30) == :ok
                    @test timedwait(() -> plot_rendered("coils"), 30) == :ok

                    raw_file = joinpath(@__DIR__, "test_files", "Koma_signal.mrd")
                    load_file!(w, :raw_data, raw_file)
                    @test !isempty(w.raw[].profiles)
                    @test timedwait(() -> w.state[] == "sig", 30) == :ok

                    # Reload reads changed samples from disk, not the previously loaded object.
                    reload_raw = joinpath(mktempdir(), "reload.mrd")
                    cp(raw_file, reload_raw)
                    load_file!(w, :raw_data, reload_raw)
                    changed_raw = deepcopy(w.raw[])
                    first(changed_raw.profiles).data .*= 2
                    save(ISMRMRDFile(reload_raw), changed_raw)
                    click_button("button_scanner_params")
                    @test timedwait(() -> w.state[] == "scanneparams", 30) == :ok
                    click_button("button_reload_raw")
                    @test timedwait(() -> first(w.raw[].profiles).data == first(changed_raw.profiles).data, 30) == :ok
                    @test timedwait(() -> plot_rendered("sig"), 30) == :ok
                end

                @testset "Close UI" begin
                    waiter = @async KomaMRI.keep_app_open(w)
                    app = w.window[].app
                    close(w.window[])
                    @test timedwait(() -> istaskdone(waiter), 30) == :ok
                    @test !isopen(w)
                    close(w)
                    @test timedwait(() -> process_exited(app.proc), 30) == :ok
                end
            finally
                close(w)
            end
        end
    end
end
