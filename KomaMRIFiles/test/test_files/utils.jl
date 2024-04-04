using KomaMRIBase, MAT, PrettyTables
TOLERANCE = 1e-6

get_theo_A_helper(r::RF) = sum(abs.(r.A)) == 0 ? [0; 0; zeros(2*length(r.A)); 0] : [0; 0; transpose([r.A r.A])[:]; 0]

# Auxiliar functions
function get_theo_t_aux(rf::RF)
    Namp, Ntim = length(rf.A), length(rf.T)
    if Namp == 1 && Ntim == 1
        return KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
    elseif Namp > 1 && Ntim == 1
        amps = KomaMRIBase.get_theo_t(rf; max_rf_samples=Inf)[2:end]
        return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
    end
    return []
end
function get_theo_A_aux(rf::RF)
    Namp, Ntim = length(rf.A), length(rf.T)
    if Namp == 1 && Ntim == 1
        return γ*get_theo_A_helper(rf)[2:end]
    elseif Namp > 1 && Ntim == 1
        amps = γ*get_theo_A_helper(rf)[2:end]
        return [amps[1]; amps[2:end-1][[i for i in 1:length(amps)-1 if i % 2 == 0]]; amps[end]]
    end
    return []
end
function get_theo_t_aux(gr::Grad)
    if !(gr.A isa Vector) && !(gr.T isa Vector)
        return KomaMRIBase.get_theo_t(gr)
    else
        return KomaMRIBase.get_theo_t(gr)[2:end]
    end
end
function get_theo_A_aux(gr::Grad)
    if !(gr.A isa Vector) && !(gr.T isa Vector)
        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)
    else
        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)[2:end]
    end
end
function get_theo_t_aux(adc::ADC)
    return (adc.N == 1) ? [adc.T/2] .+ adc.delay : [range(0, adc.T; length=adc.N)...] .+ adc.delay
end

function koma_read(filename)
    path = @__DIR__
    seq = read_seq("$(path)/pulseq_read_comparison/$(filename).seq")
    N = length(seq)
    T0 = KomaMRIBase.get_block_start_times(seq)
    return Dict(
        "rf" => Dict(
                "blocks"     => [i for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
                "times"      => [get_theo_t_aux(seq[i].RF[1]) .+ T0[i] for i in 1:N if KomaMRIBase.is_RF_on(seq[i])],
                "amplitudes" => [get_theo_A_aux(seq[i].RF[1]) for i in 1:N if KomaMRIBase.is_RF_on(seq[i])]
            ),
        "gx" => Dict(
                "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
                "times"      => [get_theo_t_aux(seq.GR[1,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])],
                "amplitudes" => [get_theo_A_aux(seq.GR[1,i]) for i in 1:N if KomaMRIBase.is_Gx_on(seq[i])]
            ),
        "gy" => Dict(
                "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
                "times"      => [get_theo_t_aux(seq.GR[2,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])],
                "amplitudes" => [get_theo_A_aux(seq.GR[2,i]) for i in 1:N if KomaMRIBase.is_Gy_on(seq[i])]
            ),
        "gz" => Dict(
                "blocks"     => [i for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
                "times"      => [get_theo_t_aux(seq.GR[3,i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])],
                "amplitudes" => [get_theo_A_aux(seq.GR[3,i]) for i in 1:N if KomaMRIBase.is_Gz_on(seq[i])]
            ),
        "adc" => Dict(
                "blocks"     => [i for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
                "times"      => [get_theo_t_aux(seq.ADC[i]) .+ T0[i] for i in 1:N if KomaMRIBase.is_ADC_on(seq[i])],
            )
    )
end

function pulseq_read(filename)
    path = @__DIR__
    seq_pulseq = matread("$(path)/pulseq_read_comparison/$(filename).mat")
    return Dict(
        "rf" => Dict(
                "blocks"     => Int.(vec(seq_pulseq["rfBlocks"])),
                "times"      => vec([vec(blockRfTimes) for blockRfTimes in seq_pulseq["rfTimes"]]),
                "amplitudes" => vec([vec(blockRfAmplitudes) for blockRfAmplitudes in seq_pulseq["rfAmplitudes"]])
            ),
        "gx" => Dict(
                "blocks"     => Int.(vec(seq_pulseq["gxBlocks"])),
                "times"      => vec([vec(blockGxTimes) for blockGxTimes in seq_pulseq["gxTimes"]]),
                "amplitudes" => vec([vec(blockGxAmplitudes) for blockGxAmplitudes in seq_pulseq["gxAmplitudes"]])
            ),
        "gy" => Dict(
                "blocks"     => Int.(vec(seq_pulseq["gyBlocks"])),
                "times"      => vec([vec(blockGyTimes) for blockGyTimes in seq_pulseq["gyTimes"]]),
                "amplitudes" => vec([vec(blockGyAmplitudes) for blockGyAmplitudes in seq_pulseq["gyAmplitudes"]])
            ),
        "gz" => Dict(
                "blocks"     => Int.(vec(seq_pulseq["gzBlocks"])),
                "times"      => vec([vec(blockGzTimes) for blockGzTimes in seq_pulseq["gzTimes"]]),
                "amplitudes" => vec([vec(blockGzAmplitudes) for blockGzAmplitudes in seq_pulseq["gzAmplitudes"]])
            ),
        "adc" => Dict(
                "blocks"     => Int.(vec(seq_pulseq["adcBlocks"])),
                "times"      => vec([vec(blockAdcTimes) for blockAdcTimes in seq_pulseq["adcTimes"]]),
            )
    )
end

function read_comparison(filename)
    # Define some names of the sequence dictionaries
    event_names = ["rf", "gx", "gy", "gz", "adc"]
    sample_names = ["times", "amplitudes"]
    # Read data for Koma and Pulseq
    seq_koma = koma_read(filename)      # it reads the .seq
    seq_pulseq = pulseq_read(filename)  # it reads the .mat (ground truth)
    # Checks for the blocks IDs of each event
    for event in event_names
        blocks_koma = seq_koma[event]["blocks"]
        blocks_pulseq = seq_pulseq[event]["blocks"]
        same_number_of_blocks = length(blocks_koma) == length(blocks_pulseq)
        @test same_number_of_blocks
        if same_number_of_blocks
            same_blocks = all(blocks_koma .== blocks_pulseq)
            @test same_blocks
            if !same_blocks
                println("Mismatch for $(event) blocks")
                data = [blocks_koma blocks_pulseq]
                pretty_table(data; header=["koma", "pulseq"])
            end
        else
            println("Length mismatch for $(event):")
            println("Koma $(length(seq_koma[event])) != Pulseq $(length(seq_pulseq[event]))")
        end
    end
    # Checks for the timings of each event
    for event in event_names
        for sample in sample_names
            if !(event == "adc" && sample == "amplitudes")
                for i in 1:length(seq_pulseq[event]["blocks"]) #1:5
                    samples_koma = seq_koma[event][sample][i]
                    samples_pulseq = seq_pulseq[event][sample][i]
                    same_points = all(.≈(abs.(samples_koma .- samples_pulseq), 0, atol=TOLERANCE))
                    @test same_points
                    if !same_points
                        println("Mismatch for $(event) $(sample). Block $(seq_pulseq[event]["blocks"][i])")
                        data = [samples_koma samples_pulseq]
                        pretty_table(data; header=["koma", "pulseq"])
                        # This is just an example about how to debug locally when there are problems
                        #if event == "gx" && sample == "amplitudes" #&& seq_pulseq[event]["blocks"][i] == 3
                        #    pretty_table(data; header=["koma", "pulseq"])
                        #    pretty_table([seq_koma[event]["times"][i] seq_pulseq[event]["times"][i]]; header=["koma_t", "pulseq_t"])
                        #    #pretty_table(data[end-10:end,:]; header=["koma", "pulseq"])
                        #end
                    end
                end
            end
        end
    end
end
