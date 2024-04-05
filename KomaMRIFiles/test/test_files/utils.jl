using KomaMRIBase, MAT, PrettyTables
TOLERANCE = 1e-6

# Auxiliar functions
#function get_theo_t_aux(rf::RF)
#    if !(rf.A isa Vector) && !(rf.T isa Vector)
#        return KomaMRIBase.get_theo_t(rf, :A)
#    else
#        return KomaMRIBase.get_theo_t(rf, :A)[2:end]
#    end
#end
#function get_theo_A_aux(rf::RF)
#    if !(rf.A isa Vector) && !(rf.T isa Vector)
#        return γ*KomaMRIBase.get_theo_A(rf, :A)
#    else
#        return γ*KomaMRIBase.get_theo_A(rf, :A)[2:end]
#    end
#end
#function get_theo_t_aux(gr::Grad)
#    if !(gr.A isa Vector) && !(gr.T isa Vector)
#        return KomaMRIBase.get_theo_t(gr)
#    else
#        return KomaMRIBase.get_theo_t(gr)[2:end]
#    end
#end
#function get_theo_A_aux(gr::Grad)
#    if !(gr.A isa Vector) && !(gr.T isa Vector)
#        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)
#    else
#        return 1e-3*γ*KomaMRIBase.get_theo_A(gr; off_val=0)[2:end]
#    end
#end
#function get_theo_t_aux(adc::ADC)
#    return (adc.N == 1) ? [adc.T/2] .+ adc.delay : [range(0, adc.T; length=adc.N)...] .+ adc.delay
#end

# For transforming Dictionary to NamedTuple recursively
namedtuple(x) = x[:]
namedtuple(d::Dict) = (; (Symbol(k) => namedtuple(v) for (k,v) in d)...)
rename_df_key(dict::Dict) = begin
    dict["Δf"] = dict["df"]
    delete!(dict, "df")
    return dict
end

function read_comparison_local(filename)
    path = @__DIR__
    event_names = [:rf, :Δf, :gx, :gy, :gz, :adc]
    type_names = [:t, :A]
    seq_koma = read_seq("$(path)/pulseq_read_comparison/$(filename).seq")
    seq_pulseq = namedtuple.(rename_df_key.(matread("$(path)/pulseq_read_comparison/$(filename).mat")["sequence"]))
    for i in 1:length(seq_koma)
        block_koma = get_samples(seq_koma, i)
        block_pulseq = seq_pulseq[i]
        for event in event_names
            for type in type_names
                samples_koma = getproperty(getproperty(block_koma, event), type)
                samples_pulseq = getproperty(getproperty(block_pulseq, event), type)
                println("Block $(i): $(event).$(type)")
                println("samples_koma: $(samples_koma)")
                println("samples_pulseq: $(samples_pulseq)")
                same_samples = all(.≈(samples_koma, samples_pulseq, atol=TOLERANCE))
                @test same_samples
                if !same_samples
                    println("Mismatch for block $(i): $(event).$(type)")
                    data = [samples_koma samples_pulseq]
                    pretty_table(data; header=["koma", "pulseq"])
                end
            end
        end
    end

end

function read_comparison(filename)
    path = @__DIR__
    seq_koma = read_seq("$(path)/pulseq_read_comparison/$(filename).seq")
    seq_pulseq = namedtuple.(rename_df_key.(matread("$(path)/pulseq_read_comparison/$(filename).mat")["sequence"]))
    for i in 1:length(seq_koma)
        block_koma = get_samples(seq_koma, i)
        block_pulseq = seq_pulseq[i]
        for (event_koma, event_pulseq) in zip(block_koma, NamedTuple{keys(block_koma)}(block_pulseq))
            @test all(.≈(event_koma.t, event_pulseq.t, atol=TOLERANCE))
            @test all(.≈(event_koma.A, event_pulseq.A, atol=TOLERANCE))
        end
    end
end
