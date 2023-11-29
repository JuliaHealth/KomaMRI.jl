using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:io in ti.tags) #verbose=true

@testitem "IO" tags=[:io] begin
    using Suppressor
    using KomaMRICore   # For RawAcquisitionData

    # Test Pulseq
    @testset "Pulseq" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        seq = @suppress read_seq(path*"/test_files/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        # Test ISMRMRD
        fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
        raw = RawAcquisitionData(fraw)
        @test raw.params["protocolName"] == "epi"
        @test raw.params["institutionName"] == "Pontificia Universidad Catolica de Chile"
        @test raw.params["encodedSize"] ≈ [101, 101, 1]
        @test raw.params["reconSize"] ≈ [102, 102, 1]
        @test raw.params["patientName"] == "brain2D_axial"
        @test raw.params["trajectory"] == "other"
        @test raw.params["systemVendor"] == "KomaMRI.jl"

        # Test signal_to_raw_data
        signal1 = Vector()
        for i=1:length(raw.profiles)
            signal1 = [signal1; raw.profiles[i].data]
        end
        rawmrd = signal_to_raw_data(signal1, seq)
        @test rawmrd.params["institutionName"] == raw.params["institutionName"]
        io = IOBuffer()
        show(io, "text/plain", rawmrd)
        @test occursin("RawAcquisitionData[", String(take!(io)))

    end
    # Test JEMRIS
    @testset "JEMRIS" begin
        path = @__DIR__
        obj = read_phantom_jemris(path*"/test_files/column1d.h5")
        @test obj.name == "column1d.h5"
    end
    # Test JEMRIS
    @testset "MRiLab" begin
        path = @__DIR__
        filename = path * "/test_files/brain_mrilab.mat"
        FRange_filename = path * "/test_files/FRange.mat" #Slab within slice thickness
        obj = read_phantom_MRiLab(filename; FRange_filename)
        @test obj.name == "brain_mrilab.mat"
    end
end
