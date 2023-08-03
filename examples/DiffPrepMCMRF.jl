# Code used to generate moment-compensated diffusion gradient waveforms
# Sequence optimization for diffusion prepared motion-compensated MRF 

using KomaMRI, KomaMRICore, JuMP, Ipopt, Dates
using LinearAlgebra: I, Bidiagonal, norm, Diagonal, Tridiagonal
using Printf
using PlotlyJS

## Aux functions
""""Calculates the normalized moments Mₖ = 1/tᵏ ∫ᵗG(τ)τᵏ dτ at the end of the sequence. """
function get_Mmatrix(seq::Sequence; axis=1, τ_sample=dur(seq))
    τ = τ_sample # Seq Duration [ms]
    T0 = cumsum([0; seq.DUR])
    M0, M1, M2, M3 = Float64[], Float64[], Float64[], Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay #Durations of delay [s]
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) #Durations of pulse [s]
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations - P0 model
            # append!(M0, δ/τ)
            # append!(M1, δ.*(T .+ δ/2)/τ^2)
            # append!(M2, δ.*(T.^2 .+ T.*δ .+ δ.^2/3)/τ^3)
            # append!(M3, δ.*(T.^3 .+ 3/2 * T.^2 .*δ .+ T.*δ.^2 .+ δ.^3/4)/τ^4)
            #Moment calculations - P1 model
            append!(M0, δ/τ^1)
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

"""Maxwell matrix."""
function get_MXmatrix(seq::Sequence; axis = 1)
    MX = Tridiagonal{Float64}[]
    τ = dur(seq) # Seq Duration [ms]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        if N > 1
            Δt = ones(N) * Gi.T / (N-1)
            dd = Δt[2:end]/6
            d = 2Δt/3
            MXi = Tridiagonal(dd, d, dd)
            push!(MX, MXi/τ^3)
        end
    end
    MX
end

"""Eddy current matrix: dG/dt * e^{-t/λ}."""
function get_ECmatrix(seq::Sequence; axis = 1, λ = 80e-3, τ_sample=dur(seq))
    T0 = cumsum([0; seq.DUR])
    EC = Float64[]
    SR = get_SRmatrix(seq)
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay #Durations of delay [s]
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) #Durations of pulse [s]
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations - P1 model
            ec = λ .* exp.((T .- τ_sample)/λ) .* (exp.(δ/λ) .- 1)
            append!(EC, - ec' * SR[i] )
        end
    end
    EC'
end

"""Eddy current matrix: dG/dt * e^{-t/λ}."""
function get_ECmatrixM0(seq::Sequence; axis = 1, λ = 80e-3, τ_sample=dur(seq))
    τ = dur(seq) # Seq Duration [ms]
    T0 = cumsum([0; seq.DUR])
    EC = Float64[]
    SR = get_SRmatrix(seq)
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay #Durations of delay [s]
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) #Durations of pulse [s]
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations - P1 model
            ec = δ * λ .+ λ^2 .* ( exp.(-(τ_sample .- T)/λ) .- exp.(-(τ_sample .- (T .+δ))/λ) ) 
            append!(EC, - ec' * SR[i])
        end
    end
    EC'
end

"Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m]."
get_Bmatrix(seq::Sequence; axis=1) = begin
    T0 = cumsum([0; seq.DUR[1:end-1]])
    #Calculating timings
    T = Float64[]
    δ = Float64[]
    for (i, s) = enumerate(seq)
        #Gradient
        Gi = s.GR[axis]
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
function write_diffprep_fwf(G1, G2, G3, bmax, Gmax, Smax; filename="./qte_vectors_input.txt", name="Maxwell2", 
    precision::Int=6, dwell_time=6.4e-6, verbose=false)
    open(filename, "w") do io
        t1 = range(0, G1.GR.dur[1] - maximum(G1.GR.delay), step=dwell_time) #length=δ2N(maximum(G1.GR.T))) #step=dwell_time) #
		t2 = range(0, G2.GR.dur[1] - maximum(G2.GR.delay), step=dwell_time) #length=δ2N(maximum(G2.GR.T)))
        t3 = range(0, G3.GR.dur[1] - maximum(G3.GR.delay), step=dwell_time) #length=δ2N(maximum(G3.GR.T)))
        maxN = max(length(t1), length(t2), length(t3))
        Gx1, Gy1, Gz1 = KomaMRICore.get_grads(G1, Array(t1).+maximum(G1.GR.delay))
		Gx2, Gy2, Gz2 = KomaMRICore.get_grads(G2, Array(t2).+maximum(G2.GR.delay))
        Gx3, Gy3, Gz3 = KomaMRICore.get_grads(G3, Array(t3).+maximum(G3.GR.delay))
        Gx1_round = round.(Gx1 ./ Gmax, digits=precision)
        Gx2_round = round.(Gx2 ./ Gmax, digits=precision)
        Gx3_round = round.(Gx3 ./ Gmax, digits=precision)
        Gy1_round = round.(Gy1 ./ Gmax, digits=precision)
        Gy2_round = round.(Gy2 ./ Gmax, digits=precision)
        Gy3_round = round.(Gy3 ./ Gmax, digits=precision)
        Gz1_round = round.(Gz1 ./ Gmax, digits=precision)
        Gz2_round = round.(Gz2 ./ Gmax, digits=precision)
        Gz3_round = round.(Gz3 ./ Gmax, digits=precision)
        if verbose
        println("Δt1=$(t1[2]-t1[1]) $(Gx1_round[1]) $(Gx1_round[end]) $(Gy1_round[1]) $(Gy1_round[end]) $(Gz1_round[1]) $(Gz1_round[end])")
        println("Δt2=$(t2[2]-t2[1]) $(Gx2_round[1]) $(Gx2_round[end]) $(Gy2_round[1]) $(Gy2_round[end]) $(Gz2_round[1]) $(Gz2_round[end])")
        println("Δt3=$(t3[2]-t3[1]) $(Gx3_round[1]) $(Gx3_round[end]) $(Gy3_round[1]) $(Gy3_round[end]) $(Gz3_round[1]) $(Gz3_round[end])")
        end
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
        if verbose
        println("SR1 = [$SRx1, $SRy1, $SRz1]")
        println("SR2 = [$SRx2, $SRy2, $SRz2]")
        println("SR3 = [$SRx3, $SRy3, $SRz3]")
        end
        println("SR = [$(max(SRx1, SRy1, SRz1)), $(max(SRx2, SRy2, SRz2)), $(max(SRx3, SRy3, SRz3))]")
        @assert (SRx1 <= Smax) && (SRx2 <= Smax) && (SRx3 <= Smax)
        @assert (SRy1 <= Smax) && (SRy2 <= Smax) && (SRy3 <= Smax)
        @assert (SRz1 <= Smax) && (SRz2 <= Smax) && (SRz3 <= Smax)
        if verbose
        println("M01 = [$(M01[1]), $(M01[2]), $(M01[3])]")
        println("M02 = [$(M02[1]), $(M02[2]), $(M02[3])]")
        println("M03 = [$(M03[1]), $(M03[2]), $(M03[3])]")
        end
        println("M0 = $(M0.*10.0^(-precision))")
        # @assert all(M0 .== 0)
        MX1 =  [sum(floor.(Int32, Gx1_round*10^precision).^2) sum(floor.(Int32, Gy1_round*10^precision).^2) sum(floor.(Int32, Gz1_round*10^precision).^2)]
        MX2 = -[sum(floor.(Int32, Gx2_round*10^precision).^2) sum(floor.(Int32, Gy2_round*10^precision).^2) sum(floor.(Int32, Gz2_round*10^precision).^2)]
        MX3 =  [sum(floor.(Int32, Gx3_round*10^precision).^2) sum(floor.(Int32, Gy3_round*10^precision).^2) sum(floor.(Int32, Gz3_round*10^precision).^2)]
        MX = MX1 .+ MX2 .+ MX3
        if verbose
        println("MX1 = [$(MX1[1]), $(MX1[2]), $(MX1[3])]")
        println("MX2 = [$(MX2[1]), $(MX2[2]), $(MX2[3])]")
        println("MX3 = [$(MX3[1]), $(MX3[2]), $(MX3[3])]")
        end
        println("MX = $(abs.(MX).*10.0^(-2precision))")
        # @assert all(MX .== 0)
        #BVAL
        t = range(0, dur(G1+G2+G3), step=dwell_time)
        Gx, Gy, Gz = KomaMRICore.get_grads(G1-G2+G3, Array(t))
        bvalx = (2π*γ)^2 * 1e-6 * sum(cumsum(Gx * dwell_time).^2 * dwell_time)
        bvaly = (2π*γ)^2 * 1e-6 * sum(cumsum(Gy * dwell_time).^2 * dwell_time) 
        bvalz = (2π*γ)^2 * 1e-6 * sum(cumsum(Gz * dwell_time).^2 * dwell_time)
        bval = round(bvalx+bvaly+bvalz, digits=3)
        println("bval_calc = [$bvalx $bvaly $bvalz] ($bval s/mm2)")
        #Header
        N1, N2, N3 = length(t1), length(t2), length(t3)
        println("N1 = $N1 N2 = $N2 N3 = $N3")
        date = "#Generated on $(now())\r\n"
        vars =  @sprintf "%s %s %s %s %s %s %s\r\n" "#Name"*" "^(length(name)-5) "N1"*" "^(length(string(N1))-1) "N2"*" "^(length(string(N2))-1) "N3"*" "^(length(string(N3))-1) "bval"*" "^(length(string(round(bmax,digits=1)))-3) "Gmax"*" "^(length(string(round(Gmax,digits=1)))-3) "Smax"
        unit =  @sprintf "%s %s %s %s\r\n" "#"*" "^(length(name)+length(string(N1))+length(string(N2))+length(string(N3))+2)  "s/mm2"*" "^(length(string(bval))-3) "mT/m"*" "^(length(string(round(Gmax,digits=1)))-3) "T/m/s"  
        line =  @sprintf "%s %i %i %i %.1f %.1f %.1f\r\n" name N1 N2 N3 bval Gmax*1e3 Smax
        write(io, date)
        write(io, vars)
        write(io, unit)
        write(io, line)
        for i = 1:maxN
            fx1, fy1, fz1 = i ≤ length(t1) ? (Gx1_round[i], Gy1_round[i], Gz1_round[i]) : (0,0,0)
            fx2, fy2, fz2 = i ≤ length(t2) ? (Gx2_round[i], Gy2_round[i], Gz2_round[i]) : (0,0,0)
            fx3, fy3, fz3 = i ≤ length(t3) ? (Gx3_round[i], Gy3_round[i], Gz3_round[i]) : (0,0,0)
            line = @sprintf "% .6f % .6f % .6f % .6f % .6f % .6f % .6f % .6f % .6f\r\n" fx1 fy1 fz1 fx2 fy2 fz2 fx3 fy3 fz3
            write(io, line)
        end
    end
end

#Params.
dwell_time = 6.4e-6
Gmax = 62e-3 # mT/m
Smaxs = [70] # mT/m/ms
axis_to_calc = ["xyz"] # ["x", "y", "z", "xyz", "xz", "yz"] 
moment_to_calc = [0, 1] #[0, 1, 2]
# 4 : HS 50ms
# 5 : HS 55ms
# 9 : BIR4 50ms
pulses_to_calc = [4] #3,4,5]
n_dwells = 4
maxwell = true #Maxwell/concomitant gradient compensation
gap_left_ms = 0; gap_left = floor(Int64, gap_left_ms * 1e-3 / (n_dwells * dwell_time))
gap_right_ms = 0; gap_right = floor(Int64, gap_right_ms * 1e-3 / (n_dwells * dwell_time))

#Eddy currents
λ_spectra = (0.1:0.5:150) * 1e-3
ec_spectraB0 = zeros(length(λ_spectra))
ec_spectraGradM0 = zeros(length(λ_spectra))
ecc_λ = [85e-3]#[0.04, 0.1, 0.24, 0.57, 1.4, 3.4, 8.3, 20, 49, 118, 288, 700] * 1e-3
ecc_α = [1.] #[51.2*ones(3); 10.24*ones(9)] / 100.
DIF = Sequence()
DIF_ref = Sequence()

for eddy = [false], gap_left_ms = [0, 1], maxwell = [true, false]
# for eddy = [false, true], gap_left_ms = [0, 1]
gap_left = floor(Int64, gap_left_ms * 1e-3 / (n_dwells * dwell_time))
for Smax = Smaxs
for pulse_type = pulses_to_calc
##############################################################################
if pulse_type == 1
    adia = "HS2"
    δ1 = 3.4560e-3
    δ2 = 7.4944e-3
    δ3 = 3.4560e-3
    Δ1 = 13.7528e-3
    Δ2 = 31.2564e-3    
elseif pulse_type == 2
    adia = "HS2"
    δ1 = 4.7040e-3
    δ2 = 9.9968e-3
    δ3 = 4.7040e-3
    Δ1 = 15.0016e-3
    Δ2 = 35.0084e-3    
elseif pulse_type == 3
    adia = "HS2"
    δ1 = 5.9520e-3
    δ2 = 12.4928e-3
    δ3 = 5.9520e-3
    Δ1 = 16.2536e-3
    Δ2 = 38.7604e-3    
elseif pulse_type == 4
    adia = "HS2"
    δ1 = 7.2064e-3
    δ2 = 14.9952e-3
    δ3 = 7.2064e-3
    Δ1 = 17.5024e-3
    Δ2 = 42.5076e-3     
elseif pulse_type == 5
    adia = "HS2"
    δ1 = 8.4544e-3
    δ2 = 17.4912e-3
    δ3 = 8.4544e-3
    Δ1 = 18.7544e-3
    Δ2 = 46.2596e-3    
elseif pulse_type == 6
    adia = "BIR4x2_3"
    δ1 = 6.9440e-3
    δ2 = 14.4704e-3
    δ3 = 6.9440e-3
    Δ1 = 10.2648e-3
    Δ2 = 27.7684e-3     
elseif pulse_type == 7
    adia = "BIR4x2_3"
    δ1 = 8.1920e-3
    δ2 = 16.9728e-3
    δ3 = 8.1920e-3
    Δ1 = 11.5136e-3
    Δ2 = 31.5204e-3    
elseif pulse_type == 8
    adia = "BIR4x2_3"
    δ1 = 9.4400e-3
    δ2 = 19.4688e-3
    δ3 = 9.4400e-3
    Δ1 = 12.7656e-3
    Δ2 = 35.2724e-3    
elseif pulse_type == 9
    adia = "BIR4x2_3"
    δ1 = 10.6944e-3
    δ2 = 21.9712e-3
    δ3 = 10.6944e-3
    Δ1 = 14.0144e-3
    Δ2 = 39.0196e-3    
elseif pulse_type == 10
    adia = "BIR4x2_3"
    δ1 = 11.9424e-3
    δ2 = 24.4672e-3
    δ3 = 11.9424e-3
    Δ1 = 15.2664e-3
    Δ2 = 42.7716e-3    
elseif pulse_type == 11
    adia = "BIR4x2_5"
    δ1 = 5.9456e-3
    δ2 = 12.4736e-3
    δ3 = 5.9456e-3
    Δ1 = 11.2632e-3
    Δ2 = 28.76668e-3     
elseif pulse_type == 12
    adia = "BIR4x2_5"
    δ1 = 7.1936e-3
    δ2 = 14.9696e-3
    δ3 = 7.1936e-3
    Δ1 = 12.5152e-3
    Δ2 = 32.5188e-3    
elseif pulse_type == 13
    adia = "BIR4x2_5"
    δ1 = 8.4480e-3
    δ2 = 17.4720e-3
    δ3 = 8.4480e-3
    Δ1 = 13.7640e-3
    Δ2 = 36.2660e-3    
elseif pulse_type == 14
    adia = "BIR4x2_5"
    δ1 = 9.6960e-3
    δ2 = 19.9744e-3
    δ3 = 9.6960e-3
    Δ1 = 15.0128e-3
    Δ2 = 40.0180e-3    
elseif pulse_type == 15               
    adia = "BIR4x2_5"
    δ1 = 10.9440e-3
    δ2 = 22.4704e-3
    δ3 = 10.9440e-3
    Δ1 = 16.2648e-3
    Δ2 = 43.7700e-3    
elseif pulse_type == 16
    adia = "BIR4x2_3"
    δ1 = 13.1968e-3
    δ2 = 26.9696e-3
    δ3 = 13.1968e-3
    Δ1 = 16.5152e-3
    Δ2 = 46.5172e-3      
elseif pulse_type == 17
    adia = "BIR4x2_5"
    δ1 = 12.1984e-3
    δ2 = 24.9728e-3
    δ3 = 12.1984e-3
    Δ1 = 17.5136e-3
    Δ2 = 47.5156e-3
elseif pulse_type == 18
    adia = "HS2"
    δ1 = 9.7088e-3
    δ2 = 19.9936e-3
    δ3 = 9.7088e-3
    Δ1 = 20.0032e-3
    Δ2 = 50.0052e-3
end
#############################################################################
δ1_new = floor(Int64, δ1 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ2_new = floor(Int64, δ2 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ3_new = floor(Int64, δ3 / dwell_time) * dwell_time # Making the waveform match the dwell time
@assert δ1_new ≈ δ1 "δ1_new = $(δ1_new*1e3) != δ1 = $(δ1*1e3)" 
@assert δ2_new ≈ δ2 "δ2_new = $(δ2_new*1e3) != δ2 = $(δ2*1e3)" 
@assert δ3_new ≈ δ3 "δ3_new = $(δ3_new*1e3) != δ3 = $(δ3*1e3)" 
rf1 = Δ1 - δ1
rf2 = Δ2 - δ2 - Δ1
# Grads - Pre-defined RF waveforms.
N1 = floor(Int64, δ1 / (n_dwells * dwell_time)) + 1; println("N1opt = $N1")
N2 = floor(Int64, δ2 / (n_dwells * dwell_time)) + 1 # δ1/N1 = δ2/N2
N3 = floor(Int64, δ3 / (n_dwells * dwell_time)) + 1

if δ3 == 0 
    N3 = 2
    δ3 = dwell_time
    rf2 = 0
    N2 = floor(Int, N1 * δ2 / δ1)  # δ1/N1 = δ2/N2
end
println("#################### pulse_type = $pulse_type ####################")
# println("δ1=$(δ1*1e3), δ2=$(δ2*1e3), δ3=$(δ3*1e3)")
global DIF =  Sequence([Grad(x -> 1e-3, δ1, N1; delay=0)])
global DIF += Sequence([Grad(x -> 1e-3, δ2, N2; delay=rf1)])
global DIF += Sequence([Grad(x -> 1e-3, δ3, N3; delay=rf2)])
Smax_discrete = Smax * 0.999
# for i=1:3
#     δ = DIF.GR[1,i].T / (length(DIF.GR[1,i].A) - 1) * 1e6
#     println("δ_$i = $δ ms, Nδ = $([δ1_new δ2_new δ3_new][i]/([N1 N2 N3][i]-1)*1e6)")
# end

#To match the samples exactly
dt = max(δ1 / (N1-1), δ2 / (N2-1), δ3 / (N3-1))
Smax_discrete = Gmax / (dt * ceil(Int64, (Gmax / Smax) / dt))
println("Smax_discrete = ", Smax_discrete)

τ = dur(DIF) # τ/Nt = Δt => Nt = τ/Δt
durT = round(Int64, round(τ*1e3)) #For the name
# Opt matrices
B =  get_Bmatrix(DIF)  #B-value
SR = get_SRmatrix(DIF) #Slew-rate matrices
MX = get_MXmatrix(DIF) #Maxwell matrices
M =  get_Mmatrix(DIF)  #Moments
EC = zeros(length(ecc_λ), N1+N2+N3)
ECM0 = zeros(length(ecc_λ), N1+N2+N3)
for (i, λ) in enumerate(ecc_λ)
    EC[i,:] = get_ECmatrix(DIF; λ) #Eddy currents B0
    ECM0[i,:] = get_ECmatrixM0(DIF; λ) #Eddy currents
end 

#EXPERIMENT
M1 =  get_Mmatrix(DIF[1],   τ_sample=dur(DIF))  #Moments
M2 =  get_Mmatrix(DIF[1:2], τ_sample=dur(DIF))  #Moments
M1 = [M1 zeros(4, N2+N3)]
M2 = [M2 zeros(4, N3)]

for k = moment_to_calc #Number of moments to null
    seq_name = eddy ? "EC$(floor(Int64, ecc_λ[1]*1e3))_" : ""
    seq_name = maxwell ? "$(seq_name)MX_MC$(k)" : "$(seq_name)MC$(k)"  #Name of the sequnce
    seq_name = adia != "" ? "$(adia)_$(seq_name)" : seq_name       #Name of the sequnce
    seq_name *= "_$durT"
    seq_name *= gap_left_ms > 0 ? "_gap" : ""
    println("#################### $seq_name ####################")
    ## Optimazation
    Mm = M[1:k+1,:]
    Mm1 = M1[1:1,:]
    Mm2 = M2[1:1,:]
    model = Model(Ipopt.Optimizer)
    # set_optimizer_attribute(model, "constr_viol_tol", 1e-16)
    # set_optimizer_attribute(model, "acceptable_tol", 1e-16)
    set_silent(model)
    @variable(model, -Gmax <= g1[1:N1] <= Gmax, start=Gmax); #max-grads
    @variable(model, -Gmax <= g2[1:N2] <= Gmax, start=Gmax); #max-grads
    @variable(model, -Gmax <= g3[1:N3] <= Gmax, start=-Gmax); #max-grads
    @objective(model, Max, [g1;-g2;g3]'*B*[g1;-g2;g3]); #b-value
    @constraint(model, moments_final,       Mm *[g1;-g2;g3] .== 0); #moments
    # @constraint(model, moments_second_rf,   Mm2*[g1;-g2;g3] .== 0); #moments
    # @constraint(model, moments_first_rf,    Mm1*[g1;-g2;g3] .== 0); #moments
    @constraint(model, slewrate, -Smax_discrete .<= [SR[1]*g1; -SR[2]*g2; SR[3]*g3] .<= Smax_discrete); #slew rate 99.9% of the Smax
    @constraint(model, ends, [g1[1]; g2[1:1+gap_right]; g3[1]; g1[N1-gap_left:N1]; g2[N2-gap_left:N2]; g3[N3]] .== 0)
    if maxwell
        @constraint(model, concomitant, g1'*MX[1]*g1 - g2'*MX[2]*g2 + g3'*MX[3]*g3 == 0); #concomitant
    end
    if eddy
        @constraint(model, eddycurrentsB0,      ecc_α'*EC*[g1; g2; g3] .== 0);      #eddy currents B0
        @constraint(model, eddycurrentsGradM0,  ecc_α'*ECM0*[g1; -g2; g3] .== 0);   #eddy currents Grads
    end
    optimize!(model)
    gx1 = value.(g1) #retrieving solution
    gx2 = value.(g2) #retrieving solution
    gx3 = value.(g3) #retrieving solution
    global gx = [gx1; -gx2; gx3]
    # Results
    bmax = objective_value(model)
    if termination_status(model) == MOI.LOCALLY_SOLVED
        println( "Solved! 😃" )
    else
        println( "NOT Solved 😢" )
    end
    ## Solution to Sequence object (for plotting)
    t = range(-1.5*rf1/2, 1.5*rf1/2, 80)
    β = 4e2 #frequency modulation param (rad/s)
    B1 = 2*13.5e-6 * sech.(β * t)
    R1 = [RF(B1, rf1, 0, δ1);;]
    R2 = [RF(B1, rf2, 0, δ2);;]
    for axis = axis_to_calc
        if     axis == "x"
            ax = 1
            global DIF =  Sequence([Grad( gx1,δ1); Grad(0,0); Grad(0,0);;],R1)
            global DIF += Sequence([Grad( gx2,δ2); Grad(0,0); Grad(0,0);;],R2)
            global DIF += Sequence([Grad( gx3,δ3); Grad(0,0); Grad(0,0);;])
            bmax = objective_value(model)
        elseif axis == "y"
            ax = 2
            global DIF =  Sequence([Grad(0,0); Grad( gx1,δ1); Grad(0,0);;],R1)
            global DIF += Sequence([Grad(0,0); Grad( gx2,δ2); Grad(0,0);;],R2)
            global DIF += Sequence([Grad(0,0); Grad( gx3,δ3); Grad(0,0);;])
            bmax = objective_value(model)
        elseif axis == "z"
            ax = 3
            global DIF =  Sequence([Grad(0,0); Grad(0,0); Grad( gx1,δ1,0);;],R1)
            global DIF += Sequence([Grad(0,0); Grad(0,0); Grad( gx2,δ2,0);;],R2)
            global DIF += Sequence([Grad(0,0); Grad(0,0); Grad( gx3,δ3,0);;])
            bmax = objective_value(model)
        elseif axis == "xz"
            ax = 1
            global DIF =  Sequence([Grad(gx1,δ1,0); Grad(0,0); Grad(gx1,δ1,0);;],R1)
            global DIF += Sequence([Grad(gx2,δ2,0); Grad(0,0); Grad(gx2,δ2,0);;],R2)
            global DIF += Sequence([Grad(gx3,δ3,0); Grad(0,0); Grad(gx3,δ3,0);;])
            bmax = 2 * objective_value(model)
        elseif axis == "yz"
            ax = 2
            global DIF =  Sequence([Grad(0,0); Grad(gx1,δ1,0); Grad(gx1,δ1,0);;],R1)
            global DIF += Sequence([Grad(0,0); Grad(gx2,δ2,0); Grad(gx2,δ2,0);;],R2)
            global DIF += Sequence([Grad(0,0); Grad(gx3,δ3,0); Grad(gx3,δ3,0);;])
            bmax = 2 * objective_value(model)
        elseif axis == "xyz"
            ax = 1
            global DIF =  Sequence([Grad(gx1,δ1,0); Grad(gx1,δ1,0); Grad(gx1,δ1,0);;],R1)
            global DIF += Sequence([Grad(gx2,δ2,0); Grad(gx2,δ2,0); Grad(gx2,δ2,0);;],R2)
            global DIF += Sequence([Grad(gx3,δ3,0); Grad(gx3,δ3,0); Grad(gx3,δ3,0);;])
            bmax = 3 * objective_value(model)
        end
        ## TO SCANNER
        path_res = "/home/ccp/DPW/G$(floor(Int,Gmax*1e3))_SR$(ceil(Int,Smax))_$axis/"
        inv = sum(DIF[1].GR[ax].A) <= 0 #if first grdient's x area is negative, invert 
        global DIF = inv ? -DIF : DIF
        # Plots
        τ = dur(DIF) * 1e3
        p1 = plot_seq(DIF; slider=false, range=[0,τ], title="$seq_name $(round(bmax, digits=2)) s/mm2")
        p2 = plot_M0(DIF; slider=false, range=[0,τ])
        p3 = plot_M1(DIF; slider=false, range=[0,τ])
        p4 = plot_M2(DIF; slider=false, range=[0,τ])
        #Eddy
        # for (i, λ) = enumerate(λ_spectra)
        #     global ec_spectraB0[i] = get_ECmatrix(DIF; λ, τ_sample=τ*1e-3) * [gx1; gx2; gx3]
        #     global ec_spectraGradM0[i] = get_ECmatrixM0(DIF; λ, τ_sample=τ*1e-3) * [gx1; gx2; gx3]
        # end
        p5 = plot_eddy_currents(DIF, ecc_λ; α=ecc_α, slider=false, range=[0,τ])
        # p6 = plot([
        #     scatter(x=λ_spectra*1e3, y=log.(abs.(ec_spectraB0)), name="EC_B0"), 
        #     scatter(x=λ_spectra*1e3, y=log.(abs.(ec_spectraGradM0)), name="EC_GradM0")
        #     ])
        # p7 = KomaMRIPlots.plot_slew_rate(DIF; slider=false, range=[0,τ])
        p = [p1; p2 p3; p5]
        # display(p)
        savefig(p, path_res*"$seq_name.svg")
        # Write
        write_diffprep_fwf(DIF[1], DIF[2], DIF[3], bmax, Gmax, Smax; 
                filename=path_res*"$seq_name.txt", name=seq_name, verbose=false)
        println( "λ0 = $(abs(round(M[1,:]'*gx/Gmax,digits=3))), λ1 = $(abs(round(M[2,:]'*gx/Gmax,digits=3))), λ2 = $(abs(round(M[3,:]'*gx/Gmax,digits=3)))" )
        println( "MX ∫g1²-∫g2²+∫g3²=$(gx1'*MX[1]*gx1 - gx2'*MX[2]*gx2 + gx3'*MX[3]*gx3)")
        println( "b-value: $(round(bmax, digits=2)) s/mm2" )
        println( "Eddy currents B0: $(ecc_α'*EC*[gx1; gx2; gx3])")
        println( "Eddy currents GradM0: $(ecc_α'*ECM0*[gx1; -gx2; gx3])")
        println( seq_name )

        #REFERENCES
        #REF M0 TRSE
        ζ = Gmax/Smax_discrete
        G_waveform = [-Gmax; -Gmax; 0; Gmax; Gmax]
        t_waveform = [(δ2-4ζ)/2; ζ; ζ; (δ2-4ζ)/2]
        println("TRSE: $((δ1-2ζ)*1e3) $(t_waveform*1e3) $((δ3-2ζ)*1e3)")
        global DIF_ref = Sequence([Grad(Gmax, δ1-2ζ, ζ); Grad(Gmax, δ1-2ζ, ζ); Grad(Gmax, δ1-2ζ, ζ);;], R1)
        global DIF_ref += Sequence([Grad(G_waveform, t_waveform, ζ); Grad(G_waveform, t_waveform, ζ); Grad(G_waveform, t_waveform, ζ);;], R2)
        global DIF_ref += Sequence([Grad(-Gmax, δ3-2ζ, ζ); Grad(-Gmax, δ3-2ζ, ζ); Grad(-Gmax, δ3-2ζ, ζ);;])
        #Plot
        τ = dur(DIF_ref) * 1e3
        p1 = plot_seq(DIF_ref; slider=false, range=[0,τ], title="TRSE_$durT")
        p2 = plot_M0(DIF_ref; slider=false, range=[0,τ])
        p3 = plot_M1(DIF_ref; slider=false, range=[0,τ])
        p4 = plot_M2(DIF_ref; slider=false, range=[0,τ])
        p5 = plot_eddy_currents(DIF_ref, ecc_λ; α=ecc_α, slider=false, range=[0,τ])
        p = [p1; p2 p3; p5]
        # display(p)
        savefig(p, path_res*"TRSE_$durT.svg")
        #Write
        write_diffprep_fwf(DIF_ref[1], DIF_ref[2], DIF_ref[3], bmax, Gmax, Smax; 
            filename=path_res*"TRSE_$durT.txt", name="TRSE_$durT", verbose=false)
    end
end
end
end
end

println("Finished! 💃")