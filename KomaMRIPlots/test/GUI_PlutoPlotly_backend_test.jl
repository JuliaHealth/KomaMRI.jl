@testitem "PlutoPlotly" tags=[:plots] begin
    using KomaMRIBase, PlutoPlotly #Testing package extension

    @testset "GUI_seq_PlutoPlotly" begin
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

        @testset "plot_seq_PlutoPlotly" begin
            #Plot sequence
            plot_seq(seq)  #Plotting the sequence
            plot_seq(seq; width=800, height=600, slider=true, show_seq_blocks=true)
            @test true          #If the previous lines fail the test will fail
        end
    end
end
