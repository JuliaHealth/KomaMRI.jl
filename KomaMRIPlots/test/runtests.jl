using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:plots in ti.tags) #verbose=true

#GUI tests
@testitem "GUI" tags=[:plots] begin
    using KomaMRICore

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

    @testset "GUI_seq" begin

        seq = PulseDesigner.EPI_example()

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
