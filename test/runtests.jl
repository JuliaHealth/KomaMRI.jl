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

@testitem "KomaUI" tags=[:koma] begin

    using Blink, Interact

    w = KomaUI(dev_tools=true; blink_show=false)

    @testset "Open UI" begin
        @test CONT_INDEX == @js w document.getElementById("content").dataset.content
    end

    @testset "PulsesGUI" begin
        @js w document.getElementById("button_pulses_seq").click()
        @test CONT_SEQUENCE == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_pulses_kspace").click()
        @test CONT_KSPACE == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_pulses_M0").click()
        @test CONT_M0 == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_pulses_M1").click()
        @test CONT_M1 == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_pulses_M2").click()
        @test CONT_M2 == @js w document.getElementById("content").dataset.content
    end

    @testset "PhantomGUI" begin
        @js w document.getElementById("button_phantom").click()
        @test CONT_PHANTOM == @js w document.getElementById("content").dataset.content
    end

    @testset "ParamsGUI" begin
        @js w document.getElementById("button_scanner").click()
        @test CONT_SCANNER_PARAMS == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_sim_params").click()
        @test CONT_SIM_PARAMS == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_rec_params").click()
        @test CONT_REC_PARAMS == @js w document.getElementById("content").dataset.content
    end

    @testset "ReconGUI" begin
        @js w document.getElementById("button_reconstruction_absI").click()
        @test CONT_ABSI == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_reconstruction_angI").click()
        @test CONT_ANGI == @js w document.getElementById("content").dataset.content
        @js w document.getElementById("button_reconstruction_absK").click()
        @test CONT_ABSK == @js w document.getElementById("content").dataset.content
    end

    @testset "SignalGUI" begin
        @js w document.getElementById("button_sig").click()
        @test CONT_SIG == @js w document.getElementById("content").dataset.content
    end

    close(w)

end

@testitem "Auxiliar Functions" tags=[:koma] begin
    using MAT
    include(joinpath(@__DIR__, "../src/ui/ExportMATFunctions.jl"))
    @testset "ExportMATFunctions" begin
        phantom = brain_phantom2D()
        sys = Scanner()
        B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
        EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
        N = 101
        FOV = 23e-2
        EPI = PulseDesigner.EPI(FOV, N, sys)
        TE = 30e-3
        d1 = TE-dur(EPI)/2-dur(EX)
        d1 = d1 > 0 ? d1 : 0
        if d1 > 0 DELAY = Delay(d1) end
        seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI #Only add delay if d1 is positive (enough space)
        seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
        path = @__DIR__
        fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
        raw = RawAcquisitionData(fraw)
        acq = AcquisitionData(raw)
        Nx, Ny = raw.params["reconSize"][1:2]
        recParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx,Ny), :densityWeighting=>true)
        aux = reconstruction(acq, recParams)
        image  = reshape(aux.data, Nx, Ny, :)
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="all")
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="sequence")
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="phantom")
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="scanner")
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="raw")
        export_2_mat(seq, phantom, sys, raw, recParams, image, pwd(); type="image")
        @test true
    end
end
