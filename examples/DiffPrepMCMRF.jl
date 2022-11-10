# Code used to generate moment-compensated diffusion gradient waveforms
# Sequence optimization for diffusion prepared motion-compensated MRF 

using KomaMRI, JuMP, Ipopt, Dates
using LinearAlgebra: I, Bidiagonal, norm
using Printf

## Aux functions
""""Calculates the normalized moments Mₖ = 1/tᵏ ∫ᵗG(τ)τᵏ dτ at the end of the sequence. """
function get_Mmatrix(seq::Sequence; axis=1)
    τ = dur(seq) * 1e3 # Seq Duration [ms]
    T0 = cumsum([0; seq.DUR])
    M0, M1, M2, M3 = Float64[], Float64[], Float64[], Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay * 1e3 #Durations of delay [ms]
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) * 1e3 #Durations of pulse [ms]
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations - P0 model
            # append!(M0, δ/τ)
            # append!(M1, δ.*(T .+ δ/2)/τ^2)
            # append!(M2, δ.*(T.^2 .+ T.*δ .+ δ.^2/3)/τ^3)
            # append!(M3, δ.*(T.^3 .+ 3/2 * T.^2 .*δ .+ T.*δ.^2 .+ δ.^3/4)/τ^4)
            #Moment calculations - P1 model
            append!(M0, δ/τ)
            append!(M1, δ.*(T)/τ^2)
            append!(M2, δ.*(T.^2 .+ δ.^2/6)/τ^3)
            append!(M3, δ.*(T.^3 .+ T .* δ.^2/2)/τ^4)
        end
    end
    [M0'; M1'; M2'; M3']
end

"""Slew rate matrix: |SR*g| ≤ Smax."""
function get_SRmatrix(seq::Sequence; axis = 1)
    SR = Bidiagonal{Float64}[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        if N > 1
            Δt = ones(N) * Gi.T / (N-1)
            dv = Δt
            ev = Δt[1:end-1]
            SRi = Bidiagonal(-1 ./ dv, 1 ./ ev, :U)
            # SRi = [SRi[1,:]' ; SRi]; SRi[1,1] = 1/Δt[1] 
            push!(SR, SRi)
        end
    end
    SR
end

"Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m]."
get_Bmatrix(seq::Sequence; axis=1) = begin
    T0 = cumsum([0; seq.DUR[1:end-1]])
    #Calculating timings
    T = Float64[]
    δ = Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay
        if N > 1
            δi = ones(N) * Gi.T / (N-1)
            Ti = [sum(δi[1:j]) for j = 1:N-1]
            Ti = T0[i] .+ delay .+ [0; Ti] #Position of pulse
            append!(δ, δi)
            append!(T, Ti)
        end
    end
    τ = dur(seq) + δ[end]
    Nsamples = length(T)
	ij = [max(i,j) for i=1:Nsamples, j=1:Nsamples]
	α = [(i==j) ? 2/3 : 1/2 for i=1:Nsamples, j=1:Nsamples]
	b = (δ' .* δ) .* (τ .- T[ij] .- α .* δ[ij])
	b_value = (2π*γ)^2*1e-6*b # Trace of B tensor
	b_value
end

## TO SCANNER (Philips)
"""Duration in [s] => samples, with dwell-time of Δt = 6.4 μs."""
δ2N(δ) = floor(Int64, δ * 156250) + 2

"""Exports diffusion preparation waveforms for their use in the scanner."""
function write_diffprep_fwf(G1, G2, G3, bmax, Gmax, Smax; filename="./qte_vectors_input.txt", name="Maxwell2", precision::Int=5, dwell_time=6.4e-6)
    open(filename, "w") do io
        t1 = range(0, maximum(G1.GR.T), step=dwell_time) #length=δ2N(maximum(G1.GR.T)))
		t2 = range(0, maximum(G2.GR.T), step=dwell_time) #length=δ2N(maximum(G2.GR.T)))
        t3 = range(0, maximum(G3.GR.T), step=dwell_time) #length=δ2N(maximum(G3.GR.T)))
        maxN = max(length(t1), length(t2), length(t3))
        Gx1, Gy1, Gz1 = KomaMRI.get_grads(G1, Array(t1).+maximum(G1.GR.delay))
		Gx2, Gy2, Gz2 = KomaMRI.get_grads(G2, Array(t2).+maximum(G2.GR.delay))
        Gx3, Gy3, Gz3 = KomaMRI.get_grads(G3, Array(t3).+maximum(G3.GR.delay))
        Gx1_round = round.(Gx1 ./ Gmax, digits=precision)
        Gx2_round = round.(Gx2 ./ Gmax, digits=precision)
        Gx3_round = round.(Gx3 ./ Gmax, digits=precision)
        Gy1_round = round.(Gy1 ./ Gmax, digits=precision)
        Gy2_round = round.(Gy2 ./ Gmax, digits=precision)
        Gy3_round = round.(Gy3 ./ Gmax, digits=precision)
        Gz1_round = round.(Gz1 ./ Gmax, digits=precision)
        Gz2_round = round.(Gz2 ./ Gmax, digits=precision)
        Gz3_round = round.(Gz3 ./ Gmax, digits=precision)
        println("Δt1=$(t1[2]-t1[1]) $(Gx1_round[1]) $(Gx1_round[end]) $(Gy1_round[1]) $(Gy1_round[end]) $(Gz1_round[1]) $(Gz1_round[end])")
        println("Δt2=$(t2[2]-t2[1]) $(Gx2_round[1]) $(Gx2_round[end]) $(Gy2_round[1]) $(Gy2_round[end]) $(Gz2_round[1]) $(Gz2_round[end])")
        println("Δt3=$(t3[2]-t3[1]) $(Gx3_round[1]) $(Gx3_round[end]) $(Gy3_round[1]) $(Gy3_round[end]) $(Gz3_round[1]) $(Gz3_round[end])")
        M01 =  [sum(floor.(Int32, Gx1_round*10^precision)) sum(floor.(Int32, Gy1_round*10^precision)) sum(floor.(Int32, Gz1_round*10^precision))]
        M02 = -[sum(floor.(Int32, Gx2_round*10^precision)) sum(floor.(Int32, Gy2_round*10^precision)) sum(floor.(Int32, Gz2_round*10^precision))]
        M03 =  [sum(floor.(Int32, Gx3_round*10^precision)) sum(floor.(Int32, Gy3_round*10^precision)) sum(floor.(Int32, Gz3_round*10^precision))]
        M0 = M01 .+ M02 .+ M03
        Gx1_diff =  maximum(abs.(Gx1_round[2:end] .- Gx1_round[1:end-1]))
        Gx2_diff =  maximum(abs.(Gx2_round[2:end] .- Gx2_round[1:end-1]))
        Gx3_diff =  maximum(abs.(Gx3_round[2:end] .- Gx3_round[1:end-1]))
        Gy1_diff =  maximum(abs.(Gy1_round[2:end] .- Gy1_round[1:end-1]))
        Gy2_diff =  maximum(abs.(Gy2_round[2:end] .- Gy2_round[1:end-1]))
        Gy3_diff =  maximum(abs.(Gy3_round[2:end] .- Gy3_round[1:end-1]))
        Gz1_diff =  maximum(abs.(Gz1_round[2:end] .- Gz1_round[1:end-1]))
        Gz2_diff =  maximum(abs.(Gz2_round[2:end] .- Gz2_round[1:end-1]))
        Gz3_diff =  maximum(abs.(Gz3_round[2:end] .- Gz3_round[1:end-1]))
        SRx1 = Gx1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRx2 = Gx2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRx3 = Gx3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy1 = Gy1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy2 = Gy2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy3 = Gy3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz1 = Gz1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz2 = Gz2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz3 = Gz3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        println("SR1 = [$SRx1, $SRy1, $SRz1]")
        println("SR2 = [$SRx2, $SRy2, $SRz2]")
        println("SR3 = [$SRx3, $SRy3, $SRz3]")
        @assert (SRx1 <= Smax) && (SRx2 <= Smax) && (SRx3 <= Smax)
        @assert (SRy1 <= Smax) && (SRy2 <= Smax) && (SRy3 <= Smax)
        @assert (SRz1 <= Smax) && (SRz2 <= Smax) && (SRz3 <= Smax)
        println("M01 = [$(M01[1]), $(M01[2]), $(M01[3])]")
        println("M02 = [$(M02[1]), $(M02[2]), $(M02[3])]")
        println("M03 = [$(M03[1]), $(M03[2]), $(M03[3])]")
        println("M0 = $(M0.*10.0^(-precision))")
        #Header
        N1, N2, N3 = length(t1), length(t2), length(t3)
        date = "#Generated on $(now())\n"
        vars =  @sprintf "%s %s %s %s %s %s %s\n" "#Name"*" "^(length(name)-5) "N1"*" "^(length(string(N1))-2) "N2"*" "^(length(string(N2))-2) "N3"*" "^(length(string(N3))-2) "bval"*" "^(length(string(round(bmax,digits=1)))-4) "Gmax"*" "^(length(string(round(Gmax,digits=1)))-3) "Smax"
        unit =  @sprintf "%s %s %s %s\n" "#"*" "^(length(name)+length(string(N1))+length(string(N2))+length(string(N3))+2)  "s/mm2"*" "^(length(string(round(bmax,digits=1)))-4) "mT/m"*" "^(length(string(round(Gmax,digits=1)))-3) "T/m/s"  
        line =  @sprintf "%s %i %i %i %.1f %.1f %.1f\n" name N1 N2 N3 bmax Gmax*1e3 Smax
        write(io, date)
        write(io, vars)
        write(io, unit)
        write(io, line)
        for i = 1:maxN
            fx1, fy1, fz1 = i ≤ length(t1) ? (Gx1_round[i], Gy1_round[i], Gz1_round[i]) : (0,0,0)
            fx2, fy2, fz2 = i ≤ length(t2) ? (Gx2_round[i], Gy2_round[i], Gz2_round[i]) : (0,0,0)
            fx3, fy3, fz3 = i ≤ length(t3) ? (Gx3_round[i], Gy3_round[i], Gz3_round[i]) : (0,0,0)
            line = @sprintf "% .5f % .5f % .5f % .5f % .5f % .5f % .5f % .5f % .5f\n" fx1 fy1 fz1 fx2 fy2 fz2 fx3 fy3 fz3
            write(io, line)
        end
    end
end

#sla# Params.
Gmax = 62e-3 # mT/m
Smax = 100   # mT/m/ms
for pulse_type = [10, 11, 12]
##############################################################################
if pulse_type == 1
    # 35ms
    #    Delta 1 --> 9.365800 || Delta 2 --> 26.865799
    # 	 delta 1 --> 7.554600 || delta 2 --> 15.688800 || delta 3 -> 7.681400
    adia = false
    Δ1, Δ2 = 9.365800e-3, 26.865799e-3
    δ1, δ2, δ3 = 7.554600e-3, 15.688800e-3, 7.681400e-3
elseif pulse_type == 2
    # 40 ms
    #     Delta 1 --> 10.615800 || Delta 2 --> 30.615799
    #     delta 1 --> 8.804600 || delta 2 --> 18.188801 || delta 3 -> 8.931400
    adia = false
    Δ1, Δ2 = 10.615800e-3, 30.615799e-3
    δ1, δ2, δ3 = 8.804600e-3, 18.188801e-3, 8.931400e-3
elseif pulse_type == 3
    # 45 ms
    #     Delta 1 --> 11.865800 || Delta 2 --> 34.365799
    # 	  delta 1 --> 10.054600 || delta 2 --> 20.688801 || delta 3 -> 10.181400
    adia = false
    Δ1, Δ2 = 11.865800e-3, 34.365799e-3
    δ1, δ2, δ3 = 10.054600e-3, 20.688801e-3, 10.181400e-3
elseif pulse_type == 4
    # 50 ms
    #     Delta 1 --> 13.115800 || Delta 2 --> 38.115799
    #     delta 1 --> 11.304600 || delta 2 --> 23.188801 || delta 3 -> 11.431400
    adia = false
    Δ1, Δ2 = 13.115800e-3, 38.115799e-3
    δ1, δ2, δ3 = 11.304600e-3, 23.188801e-3, 11.431400e-3
elseif pulse_type == 5
    # 55 ms
    #     Delta 1 --> 14.368100 || Delta 2 --> 41.869000
    #     delta 1 --> 12.556900 || delta 2 --> 25.689700 || delta 3 -> 12.684900
    adia = false
    Δ1, Δ2 = 14.368100e-3, 41.869000e-3
    δ1, δ2, δ3 = 12.556900e-3, 25.689700e-3, 12.684900e-3
elseif pulse_type == 6
    # 60 ms
    #     Delta 1 --> 15.615800 || Delta 2 --> 45.615799
    #     delta 1 --> 13.804600 || delta 2 --> 28.188801 || delta 3 -> 13.931400
    adia = false
    Δ1, Δ2 = 15.615800e-3, 45.615799e-3
    δ1, δ2, δ3 = 13.804600e-3, 28.188801e-3, 13.931400e-3
elseif pulse_type == 7
    # 35s adiab 750deg
    #     Delta 1 --> 13.651300 || Delta 2 --> 31.155400
    #     delta 1 --> 3.270500 || delta 2 --> 7.123300 || delta 3 -> 3.398500
    adia = true
    Δ1, Δ2 = 13.651300e-3, 31.155400e-3
    δ1, δ2, δ3 = 3.270500e-3, 7.123300e-3, 3.398500e-3
elseif pulse_type == 8
    # 40s adiab 750deg
    #     Delta 1 --> 14.905700 || Delta 2 --> 34.905700
    #     delta 1 --> 4.524900 || delta 2 --> 9.619200 || delta 3 -> 4.646600
    adia = true
    Δ1, Δ2 = 14.905700e-3, 34.905700e-3
    δ1, δ2, δ3 = 4.524900e-3, 9.619200e-3, 4.646600e-3
elseif pulse_type == 9
    # 45s adiab 750deg
    #     Delta 1 --> 16.153700 || Delta 2 --> 38.656200
    #     delta 1 --> 5.772900 || delta 2 --> 12.121700 || delta 3 -> 5.900900 
    adia = true
    Δ1, Δ2 = 16.153700e-3, 38.656200e-3
    δ1, δ2, δ3 = 5.772900e-3, 12.121700e-3, 5.900900e-3
elseif pulse_type == 10
    # 50s adiab 750deg
    #     Delta 1 --> 17.401700 || Delta 2 --> 42.406600
    #     delta 1 --> 7.020900 || delta 2 --> 14.624100 || delta 3 -> 7.148900
    adia = true
    Δ1, Δ2 = 17.401700e-3, 42.406600e-3
    δ1, δ2, δ3 = 7.020900e-3, 14.624100e-3, 7.148900e-3
elseif pulse_type == 11
    # 55s adiab 750deg 
    #     Delta 1 --> 18.656100 || Delta 2 --> 46.157000
    #     delta 1 --> 8.275300 || delta 2 --> 17.120100 || delta 3 -> 8.396900
    adia = true
    Δ1, Δ2 = 18.656100e-3, 46.157000e-3
    δ1, δ2, δ3 = 8.275300e-3, 17.120100e-3, 8.396900e-3
elseif pulse_type == 12
    # 60s adiab 750deg
    #     Delta 1 --> 19.904100 || Delta 2 --> 49.907400
    #     delta 1 --> 9.523300 || delta 2 --> 19.622500 || delta 3 -> 9.651300
    adia = true
    Δ1, Δ2 = 19.904100e-3, 49.907400e-3
    δ1, δ2, δ3 = 9.523300e-3, 19.622500e-3, 9.651300e-3
end
##############################################################################

N1 = 400 # You can solve the opt problem in a lower time resolution or use δ2N(dur_grad) 
path_file = "/home/ccp/"
maxwell = true #maxwell or concomitant gradient compensation
sym = false
# Timings
dwell_time = 6.4e-6
δ1 = floor( δ1 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ2 = floor( δ2 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ3 = floor( δ3 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ3 = sym ? δ1 : δ3
rf1 = Δ1 - δ1
rf2 = Δ2 - δ2 - Δ1
# Grads - Pre-defined RF waveforms.
τ = Δ2 + δ3 # τ/Nt = Δt => Nt = τ/Δt
durT = ceil(Int64, τ*1e3) #For the name
N2 = floor(Int, (N1 - 1) * δ2 / δ1) # δ1/N1 = δ2/N2
N3 = floor(Int, (N1 - 1) * δ3 / δ1)
DIF =  Sequence([Grad(x -> 1e-3, δ1, N1; delay=0)])
DIF += Sequence([Grad(x -> 1e-3, δ2, N2; delay=rf1)])
DIF += Sequence([Grad(x -> 1e-3, δ3, N3; delay=rf2)])
# Opt matrices
B =  get_Bmatrix(DIF)  #B-value
SR = get_SRmatrix(DIF) #Slew-rate matrices
M =  get_Mmatrix(DIF)  #Moments

for k = [0, 1] #Number of moments to null
    seq_name = maxwell ? "MX_MC$(k)_$durT" : "MC$(k)_$durT"  #Name of the sequnce
    seq_name = adia ? "$(seq_name)_adia" : seq_name       #Name of the sequnce
    ## Optimazation
    Mm = M[1:k+1,:]
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, -Gmax <= g1[1:N1] <= Gmax, start=Gmax); #max-grads
    @variable(model, -Gmax <= g2[1:N2] <= Gmax, start=Gmax); #max-grads
    @variable(model, -Gmax <= g3[1:N3] <= Gmax, start=Gmax); #max-grads
    @objective(model, Max, [g1;g2;g3]'*B*[g1;g2;g3]); #b-value
    @constraint(model, moments, Mm*[g1;g2;g3] .== 0); #moments
    @constraint(model, slewrate, -Smax*0.999 .<= [SR[1]*g1; SR[2]*g2; SR[3]*g3] .<= Smax*0.999); #slew rate 99.9% of the SR
    @constraint(model, ends, [g1[1]; g2[1]; g3[1]; g1[N1]; g2[N2]; g3[N3]] .== 0)
    if maxwell
        @constraint(model, concomitant, sum(g1.^2) - sum(g2.^2) + sum(g3.^2) == 0); #concomitant
    end
    optimize!(model);
    gx1 = value.(g1) #retrieving solution
    gx2 = value.(g2) #retrieving solution
    gx3 = value.(g3) #retrieving solution
    # Results
    gx = [gx1; gx2; gx3]
    bmax = gx'*B*gx
    println( "λ0 = $(abs(round(M[1,:]'*gx/Gmax,digits=3))), λ1 = $(abs(round(M[2,:]'*gx/Gmax,digits=3))), λ2 = $(abs(round(M[3,:]'*gx/Gmax,digits=3)))" )
    println( "b-value: $(round(bmax, digits=2)) s/mm2" )
    println( seq_name )
    if termination_status(model) == MOI.LOCALLY_SOLVED
        println( "Solved! 😃"  )
    else
        println( "NOT Solved 😢"  )
    end
    ## Solution to Sequence object (for plotting)
    B1 = 15e-6
    β = 1100
    R1 = [RF(t->B1*sech(β*(t-rf1/2)), rf1; delay=δ1);;]
    R2 = [RF(t->B1*sech(β*(t-rf1/2)), rf2; delay=δ2);;]
    for axis = ["x", "y", "z"]
        if     axis == "x"
            ax = 1
            DIF =  Sequence([Grad( gx1,δ1); Grad(0,0); Grad(0,0);;],R1)
            DIF += Sequence([Grad(-gx2,δ2); Grad(0,0); Grad(0,0);;],R2)
            DIF += Sequence([Grad( gx3,δ3); Grad(0,0); Grad(0,0);;])
        elseif axis == "y"
            ax = 2
            DIF =  Sequence([Grad(0,0); Grad( gx1,δ1); Grad(0,0);;],R1)
            DIF += Sequence([Grad(0,0); Grad(-gx2,δ2); Grad(0,0);;],R1)
            DIF += Sequence([Grad(0,0); Grad( gx3,δ3); Grad(0,0);;])
        elseif axis == "z"
            ax = 3
            DIF =  Sequence([Grad(0,0); Grad(0,0); Grad( gx1,δ1,0);;],R1)
            DIF += Sequence([Grad(0,0); Grad(0,0); Grad(-gx2,δ2,0);;],R2)
            DIF += Sequence([Grad(0,0); Grad(0,0); Grad( gx3,δ3,0);;])
        end
        ## TO SCANNER
        path_res = "/home/ccp/DiffPrepWaveforms/G$(floor(Int,Gmax*1e3))_SR$(floor(Int,Smax))_$axis/"
        inv = DIF[1].GR[ax].A[2] <= 0 #if first grdient's x component goes down, invert 
        DIFinv = inv ? -DIF : DIF
        # Write
        write_diffprep_fwf(DIFinv[1], DIFinv[2], DIFinv[3], bmax, Gmax, Smax; filename=path_res*"$seq_name.txt", name=seq_name)
        # Plots
        R90 = RF(B1, 0.35e-3)
        p = plot_seq(R90+DIFinv+R90; darkmode=false, slider=false, range=[-1 dur(DIFinv)*1e3+1])
        savefig(p, path_res*"$seq_name.svg")
    end
end
end
println("Finished! 💃")