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

#"""
#Execute function f() with a timeout of `timeout` seconds. Returns the
#result of f() or `nothing` in the case of a timeout.
#"""
#function with_timeout(f::Function, timeout)
#    c = Channel{Any}(1)
#    @async begin
#        put!(c, f())
#    end
#    @async begin
#        sleep(timeout)
#        put!(c, nothing)
#    end
#    take!(c)
#end

@testitem "KomaUI" tags=[:koma] begin

    using Blink, Interact

#    @testset "button click" begin
#        scope = Scope()
#        obs = Observable(scope, "obs", false)
#        obschannel = Channel(1)
#        on((x) -> push!(obschannel, x), obs)
#        scope(dom"button#mybutton"(
#            events=Dict(
#                "click" => @js function()
#                    $obs[] = true
#                end
#            )
#        ))
#        w = Window(Dict(:show => false); body=scope)
#
#        # Sleep to allow WebIO scope to mount in Electron
#        sleep(0.25)
#
#        @js w document.querySelector("#mybutton").click()
#        did_click = with_timeout(() -> take!(obschannel), 5)
#        @test did_click
#    end

#    include(joinpath(@__DIR__, "../src/reconstruction/Recon.jl"))
#
#    w = Window(Blink.Dict(:show => false), async=false);
#    body!(w, """<div class="container" style="padding: 0px !important;" id="content">""");
#
#    phantom = brain_phantom2D()
#    sys = Scanner()
#    B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
#    EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
#    N = 101
#    FOV = 23e-2
#    EPI = PulseDesigner.EPI(FOV, N, sys)
#    TE = 30e-3
#    d1 = TE-dur(EPI)/2-dur(EX)
#    d1 = d1 > 0 ? d1 : 0
#    if d1 > 0 DELAY = Delay(d1) end
#    seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI #Only add delay if d1 is positive (enough space)
#    seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
#    simParams = KomaMRICore.default_sim_params()
#
#    path = @__DIR__
#    fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
#    darkmode = true
#    raw_ismrmrd = RawAcquisitionData(fraw)
#
#    acq = AcquisitionData(raw_ismrmrd)
#    Nx, Ny = raw_ismrmrd.params["reconSize"][1:2]
#    recParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx,Ny), :densityWeighting=>true)
#    aux = reconstruction(acq, recParams)
#    image  = reshape(aux.data, Nx, Ny, :)
#    kspace = KomaMRI.fftc(image)
#
#    seq_obs = Observable{Sequence}(seq)
#    pha_obs = Observable{Phantom}(phantom)
#    sig_obs = Observable{RawAcquisitionData}(raw_ismrmrd)
#    img_obs = Observable{Any}(image)
#
    @testset "PulsesGUI" begin
        w = KomaUI(dev_tools=true)
        sleep(0.25)
        @js w document.getElementById("button_pulses_seq").click()
        @js w document.getElementById("button_pulses_kspace").click()
        @js w document.getElementById("button_pulses_M0").click()
        @js w document.getElementById("button_pulses_M1").click()
        @js w document.getElementById("button_pulses_M2").click()
        close(w)
        @test true
        #try
        #    w = KomaUI(dev_tools=true)
        #    @js w Blink.msg("pulses_seq", 1)
        #    @js w Blink.msg("pulses_kspace", 1)
        #    @js w Blink.msg("pulses_M0", 1)
        #    @js w Blink.msg("pulses_M1", 1)
        #    @js w Blink.msg("pulses_M2", 1)
        #    close(w)
        #    @test true
        #catch err
        #    include(joinpath(@__DIR__, "../src/ui/PulsesGUI_seq.jl"))
        #    include(joinpath(@__DIR__, "../src/ui/PulsesGUI_kspace.jl"))
        #    include(joinpath(@__DIR__, "../src/ui/PulsesGUI_M0.jl"))
        #    include(joinpath(@__DIR__, "../src/ui/PulsesGUI_M1.jl"))
        #    include(joinpath(@__DIR__, "../src/ui/PulsesGUI_M2.jl"))
        #    @test true
        #end
    end

    @testset "PhantomGUI" begin
        w = KomaUI(dev_tools=true)
        sleep(0.25)
        @js w document.getElementById("button_phantom").click()
        close(w)
        @test true
#        #try
#        #    w = KomaUI(dev_tools=true)
#        #    @js w Blink.msg("phantom", 1)
#        #    close(w)
#        #    @test true
#        #catch err
#            include(joinpath(@__DIR__, "../src/ui/PhantomGUI.jl"))
#            @test true
#        #end
    end

    @testset "ParamsGUI" begin
        w = KomaUI(dev_tools=true)
        sleep(0.25)
        @js w document.getElementById("button_scanner").click()
        @js w document.getElementById("button_sim_params").click()
        @js w document.getElementById("button_rec_params").click()
        close(w)
        @test true
#        #try
#        #    w = KomaUI(dev_tools=true)
#        #    @js w Blink.msg("scanner", 1)
#        #    @js w Blink.msg("sim_params", 1)
#        #    @js w Blink.msg("rec_params", 1)
#        #    close(w)
#        #    @test true
#        #catch err
#            include(joinpath(@__DIR__, "../src/ui/ScannerParams_view.jl"))
#            include(joinpath(@__DIR__, "../src/ui/SimParams_view.jl"))
#            include(joinpath(@__DIR__, "../src/ui/RecParams_view.jl"))
#            @test true
#        #end
#
    end

    @testset "ReconGUI" begin
        w = KomaUI(dev_tools=true)
        sleep(0.25)
        @js w document.getElementById("button_reconstruction_absI").click()
        @js w document.getElementById("button_reconstruction_angI").click()
        @js w document.getElementById("button_reconstruction_absK").click()
        close(w)
        @test true
#        #try
#        #    w = KomaUI(dev_tools=true)
#        #    @js w Blink.msg("reconstruction_absI", 1)
#        #    @js w Blink.msg("reconstruction_angI", 1)
#        #    @js w Blink.msg("reconstruction_absK", 1)
#        #    close(w)
#        #    @test true
#        #catch err
#            include(joinpath(@__DIR__, "../src/ui/ReconGUI_absI.jl"))
#            include(joinpath(@__DIR__, "../src/ui/ReconGUI_angI.jl"))
#            include(joinpath(@__DIR__, "../src/ui/ReconGUI_absK.jl"))
#            @test true
#        #end
    end

    @testset "SignalGUI" begin
        w = KomaUI(dev_tools=true)
        sleep(0.25)
        @js w document.getElementById("button_sig").click()
        close(w)
        @test true
#        #try
#        #    w = KomaUI(dev_tools=true)
#        #    @js w Blink.msg("sig", 1)
#        #    close(w)
#        #    @test true
#        #catch err
#            include(joinpath(@__DIR__, "../src/ui/SignalGUI.jl"))
#            @test true
#        #end
    end

    @testset "Open UI" begin
        KomaUI()
        sleep(0.25)
        @test true
    end
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
