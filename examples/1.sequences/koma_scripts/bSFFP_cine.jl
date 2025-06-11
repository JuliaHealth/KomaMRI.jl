include("bSSFP.jl")

function bSSFP_cine(
    FOV::Float64, 
	N_matrix::Int, 
	TR::Float64, 
	flip_angle, 
	RRs,
	N_phases,
	sys::Scanner; 
	G=[0,0,0], 
	Δf=0,
    adc_duration=1e-3,
    N_dummy_cycles=10
)
    function RR(i)
        idx = (i - 1) % length(RRs) + 1
        return RRs[idx]
    end

    seq = Sequence()
    base_seq =  bSSFP(FOV, N_matrix, TR, flip_angle, sys; Δf=Δf, adc_duration=adc_duration)

    n = 1
    for i in 0:(N_matrix - 1) # 1 vps (Views Per Segment)
        line = base_seq[6*i .+ (1:6)]

        if N_dummy_cycles > 0 && i == 0
            dummy_dur = N_dummy_cycles * dur(line)
            rr_sum = RR(n)
            while dummy_dur > rr_sum
                n += 1
                rr_sum += RR(n)
                
            end
            d_pre = rr_sum - dummy_dur
            seq += Delay(d_pre)
        end
        
        for j in 1:(N_phases + N_dummy_cycles)
            l = copy(line)

            l[1].RF[1].A  *= (-1)^(j)         # Sign of the RF pulse is alteranted every TR
            l[4].ADC[1].ϕ  = j%2==0 ? 0 : π   # so, phase of the ADC is consecuently alterned between 0 and π

            if j in (1:N_dummy_cycles)
                l = 0*l
                l[4].ADC[1].N = 0
            end

            seq += l
        end

        phases_and_dummy_dur = (N_phases + N_dummy_cycles)* dur(line)
        n += 1
        rr_sum = RR(n)
        while phases_and_dummy_dur > rr_sum
            n += 1
            rr_sum += RR(n)
           
        end
        d_post = rr_sum - phases_and_dummy_dur
        seq += Delay(d_post)
    end

    return seq
end