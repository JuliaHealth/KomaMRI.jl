#GUI tests
@testitem "PlotlyJS" tags=[:plots] begin
    using KomaMRIBase, MRIFiles

    @testset "GUI_phantom" begin
        ph = brain_phantom2D()    #2D phantom

        @testset "plot_phantom_map_rho" begin
            plot_phantom_map(ph, :ρ, width=800, height=600) #Plotting the phantom's rho map
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

        @testset "plot_phantom_map_x" begin
            plot_phantom_map(ph, :x) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_phantom_map_w" begin
            plot_phantom_map(ph, :Δw) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_phantom_map_2dview" begin
            plot_phantom_map(ph, :ρ, view_2d=true) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end
    end

    @testset "GUI_motion_phantom" begin
        ph = brain_phantom2D()    #2D phantom
        ph.motion = MotionList(Translate(0.1, 0.1, 0.1, TimeRange(0.0, 1.0), SpinRange(1:1000)))

        @testset "plot_motion_phantom_map_rho" begin
            plot_phantom_map(ph, :ρ, width=800, height=600, max_spins=1_000) #Plotting the phantom's rho map (set max_spins=1_000)
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_motion_phantom_map_T1" begin
            plot_phantom_map(ph, :T1) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_motion_phantom_map_T2" begin
            plot_phantom_map(ph, :T2) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_motion_phantom_map_x" begin
            plot_phantom_map(ph, :x) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_motion_phantom_map_w" begin
            plot_phantom_map(ph, :Δw) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end

        @testset "plot_motion_phantom_map_2dview" begin
            plot_phantom_map(ph, :ρ, view_2d=true) #Plotting the phantom's rho map
            @test true                #If the previous line fails the test will fail
        end
    end

    @testset "GUI_seq" begin
        #KomaCore definition of a sequence:
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
            plot_seq(seq; width=800, height=600, slider=true, show_seq_blocks=true)
            @test true          #If the previous lines fail the test will fail
        end

        @testset "plot_kspace" begin
            #Plot k-space
            plot_kspace(seq; width=800, height=600) #Plotting the k-space
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_M0" begin
            #Plot M0
            plot_M0(seq)        #Plotting the M0
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_M1" begin
            #Plot M1
            plot_M1(seq)        #Plotting the M0
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_M2" begin
            #Plot M2
            plot_M2(seq)        #Plotting the M2
            @test true          #If the previous line fails the test will fail
        end

        @testset "plot_eddy_currents" begin
            #Plot M2
            plot_eddy_currents(seq, 80e-3) #Plotting the plot_eddy_currents
            @test true              #If the previous line fails the test will fail
        end

        @testset "plot_slew_rate" begin
            plot_slew_rate(seq)
            @test true
        end

        @testset "plot_seqd" begin
            plot_seqd(seq)
            @test true
        end
    end

    @testset "GUI_dict_html" begin
        #Define a dictionary and Plot the dictionary table
        sys = Scanner()
        sys_dict = Dict("B0" => sys.B0,
                "B1" => sys.B1,
                "Gmax" => sys.Gmax,
                "Smax" => sys.Smax,
                "ADC_dt" => sys.ADC_Δt,
                "seq_dt" => sys.seq_Δt,
                "GR_dt" => sys.GR_Δt,
                "RF_dt" => sys.RF_Δt,
                "RF_ring_down_T" => sys.RF_ring_down_T,
                "RF_dead_time_T" => sys.RF_dead_time_T,
                "ADC_dead_time_T" => sys.ADC_dead_time_T)
        plot_dict(sys_dict)
        @test true
    end

    @testset "GUI_signal" begin
        path = @__DIR__
        fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
        raw = RawAcquisitionData(fraw)
        plot_signal(raw, width=800, height=600)
        @test true          #If the previous line fails the test will fail
    end

    @testset "GUI_recon" begin
        #???
    end

end
