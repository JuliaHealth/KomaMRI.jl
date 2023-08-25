function export_2_mat_sequence(seq, matfolder; matfilename="seq_sequence.mat")
	max_rf_samples=100
    N = length(seq)
    ΔT = KomaMRICore.durs(seq)
    T0 = cumsum([0; ΔT],dims=1)
    off_val = Inf #This removes the unnecessary points in the plot

    #GRADS
    t1x = vcat([KomaMRICore.get_theo_t(seq.GR[1,i]) .+ T0[i] for i=1:N]...)
    t1y = vcat([KomaMRICore.get_theo_t(seq.GR[2,i]) .+ T0[i] for i=1:N]...)
    t1z = vcat([KomaMRICore.get_theo_t(seq.GR[3,i]) .+ T0[i] for i=1:N]...)
    Gx =  vcat([KomaMRICore.get_theo_A(seq.GR[1,i]) for i=1:N]...)
    Gy =  vcat([KomaMRICore.get_theo_A(seq.GR[2,i]) for i=1:N]...)
    Gz =  vcat([KomaMRICore.get_theo_A(seq.GR[3,i]) for i=1:N]...)
    GRADS = hcat(t1x, t1y, t1z, Gx, Gy, Gz)
    #RFS
    t2 =  vcat([KomaMRICore.get_theo_t(seq.RF[1,i];max_rf_samples) .+ T0[i] for i=1:N]...)
    R =   vcat([KomaMRICore.get_theo_A(r;off_val,max_rf_samples) for r = seq.RF]...)
    RFS = hcat(t2, R)
    #ADC
    t3 =  vcat([KomaMRICore.get_theo_t(seq.ADC[i])  .+ T0[i] for i=1:N]...)
    D =   vcat([KomaMRICore.get_theo_A(d;off_val) for d = seq.ADC]...)
    ADCS = hcat(t3, D)

    seq_dict = Dict("GRAD" => GRADS,
                    "RF" => RFS,
                    "ADC" => ADCS,
                    "DUR" => seq.DUR,
                    "DEF" => seq.DEF)
    matwrite(joinpath(matfolder, matfilename), Dict("sequence" => seq_dict))
end

function export_2_mat_kspace(seq, matfolder; matfilename="seq_kspace.mat")
    kspace, kspace_adc = get_kspace(seq; Δt=1)
    matwrite(joinpath(matfolder, matfilename), Dict("kspace" => kspace, "kspace_adc" => kspace_adc))
end

function export_2_mat_moments(seq, matfolder; matfilename="seq_moments.mat")
    dt = 1
    t, Δt = KomaMRICore.get_uniform_times(seq, dt)
    t = t[1:end-1]
    k0, _ =  KomaMRICore.get_kspace(seq; Δt=dt)
    k1, _ =  KomaMRICore.get_M1(seq; Δt=dt)
    k2, _ =  KomaMRICore.get_M2(seq; Δt=dt)
    moments = hcat(t, k0, k1, k2)
    matwrite(joinpath(matfolder, matfilename), Dict("moments" => moments))
end

function export_2_mat_phantom(phantom, matfolder; matfilename="phantom.mat")
    phantom_dict = Dict("name" => phantom.name,
                "columns" => ["x", "y", "z", "rho", "T1", "T2", "T2s", "delta_omega"],
                "data" => hcat(phantom.x, phantom.y, phantom.z, phantom.ρ, phantom.T1, phantom.T2, phantom.T2s, phantom.Δw))
    matwrite(joinpath(matfolder, matfilename), Dict("phantom" => phantom_dict))
end

function export_2_mat_scanner(sys, matfolder; matfilename="scanner.mat")
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
    matwrite(joinpath(matfolder, matfilename), Dict("scanner" => sys_dict))
end

function export_2_mat_raw(raw_ismrmrd, matfolder; matfilename="raw.mat");
    if haskey(raw_ismrmrd.params, "userParameters")
        dict_for_mat = Dict()
        dict_user_params = raw_ismrmrd.params["userParameters"]
        for (key, value) in dict_user_params
            if key == "Δt_rf"
                dict_for_mat["dt_rf"] = value
            elseif key == "Δt"
                dict_for_mat["dt_rf"] = value
            else
                dict_for_mat[key] = value
            end
        end
        matwrite(joinpath(matfolder, "sim_params.mat"), dict_for_mat)

        not_Koma = raw_ismrmrd.params["systemVendor"] != "KomaMRI.jl"
        t = Float64[]
        signal = ComplexF64[]
        current_t0 = 0
        for p in raw_ismrmrd.profiles
        	dt = p.head.sample_time_us != 0 ? p.head.sample_time_us * 1e-3 : 1
        	t0 = p.head.acquisition_time_stamp * 1e-3 #This parameter is used in Koma to store the time offset
            N =  p.head.number_of_samples != 0 ? p.head.number_of_samples : 1
            if not_Koma
        		t0 = current_t0 * dt
                current_t0 += N
            end
            if N != 1
                append!(t, t0.+(0:dt:dt*(N-1)))
            else
                append!(t, t0)
            end
            append!(signal, p.data[:,1]) #Just one coil
            #To generate gap
            append!(t, t[end])
            append!(signal, [Inf + Inf*1im])
        end
        raw_dict = hcat(t, signal)
        matwrite(joinpath(matfolder, matfilename), Dict("raw" => raw_dict))
    end

end

function export_2_mat_image(image, rec_params, matfolder; matfilename="image.mat")
    if haskey(rec_params, :reconSize)
        dict_rec_params = Dict("reco" => rec_params[:reco],
                            "Nx" => rec_params[:reconSize][1],
                            "Ny" => rec_params[:reconSize][2])
        matwrite(joinpath(matfolder, "rec_params.mat"), Dict("rec_params" => dict_rec_params))
    end

    matwrite(joinpath(matfolder, matfilename), Dict("image" => image))
end

function export_2_mat(seq, phantom, sys, raw_ismrmrd, rec_params, image, matfolder; type="all", matfilename="data.mat")
    head = splitext(matfilename)[1]
    if type=="all"
        export_2_mat_sequence(seq, matfolder)
        export_2_mat_kspace(seq, matfolder)
        export_2_mat_moments(seq, matfolder)
        export_2_mat_phantom(phantom, matfolder)
        export_2_mat_scanner(sys, matfolder)
        export_2_mat_raw(raw_ismrmrd, matfolder)
        export_2_mat_image(image, rec_params, matfolder)
    elseif type=="sequence"
		export_2_mat_sequence(seq, matfolder; matfilename=(head*"_sequence.mat"))
        export_2_mat_kspace(seq, matfolder; matfilename=(head*"_kspace.mat"))
        export_2_mat_moments(seq, matfolder; matfilename=(head*"_moments.mat"))
	elseif type=="phantom"
		export_2_mat_phantom(phantom, matfolder; matfilename=(head*"_phantom.mat"))
    elseif type=="scanner"
		export_2_mat_scanner(sys, matfolder; matfilename=(head*"_scanner.mat"))
    elseif type=="raw"
		export_2_mat_raw(raw_ismrmrd, matfolder; matfilename=(head*"_raw.mat"))
    elseif type=="image"
		export_2_mat_image(image, rec_params, matfolder; matfilename=(head*"_image.mat"))
	end

    strToast = "<ul><li><b>Name:</b> " * matfilename * "</li><li><b>Path:</b> " * matfolder * "</li></ul>"
    if type=="all"
        strToast = "<ul><li><b>Path:</b> " * matfolder * "</li></ul>"
    end
    return strToast
end
