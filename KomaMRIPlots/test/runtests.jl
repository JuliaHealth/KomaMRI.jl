using TestItems, TestItemRunner

#GUI tests
@testitem "GUI" begin
    using Suppressor, KomaMRICore
    @testset "GUI_phantom" begin
        ph = brain_phantom2D()    #2D phantom

        @testset "plot_phantom_map_rho" begin
            plot_phantom_map(ph, :ρ) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_phantom_map_T1" begin
            plot_phantom_map(ph, :T1) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_phantom_map_T2" begin
            plot_phantom_map(ph, :T2) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end
    end

    @testset "GUI_seq" begin
        #RF construction
        sys = Scanner()
        B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
        EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
        #ACQ construction
        N = 101
        FOV = 23e-2
        EPI = PulseDesigner.EPI(FOV, N, sys)
        TE = 30e-3
        d1 = TE-dur(EPI)/2-dur(EX)
        d1 = d1 > 0 ? d1 : 0
        if d1 > 0 DELAY = Delay(d1) end
        #Sequence construction
        seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI #Only add delay if d1 is positive (enough space)
        seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3

        @testset "plot_seq" begin
            #Plot sequence
            plot_seq(seq)  #Plotting the sequence
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_kspace" begin
            #Plot k-space
            plot_kspace(seq)    #Plotting the k-space
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_M0" begin
            #Plot M0
            plot_M0(seq)        #Plotting the M0
            @test true          #If the previous line fails the test will fail
        end
    end

    @testset "GUI_signal" begin
        path = @__DIR__
        fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
        raw = RawAcquisitionData(fraw)
        plot_signal(raw)
        @test true          #If the previous line fails the test will fail
    end

    @testset "GUI_recon" begin
        #???
    end

end

@run_package_tests filter=ti->!(:skipci in ti.tags)
