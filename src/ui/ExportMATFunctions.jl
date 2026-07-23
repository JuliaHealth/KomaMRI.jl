function export_2_mat_sequence(seq, matfolder; matfilename="seq_sequence.mat")

    # Get the sequence event waveforms
    samples = get_samples(seq)

    # Create a dictionary with the sampled event waveforms
    Gx = hcat(samples.gx.t, samples.gx.A)
    Gy = hcat(samples.gy.t, samples.gy.A)
    Gz = hcat(samples.gz.t, samples.gz.A)
    RF_AM = hcat(samples.rf.t, samples.rf.A)
    RF_FM = hcat(samples.Δf.t, samples.Δf.A)
    ADCS = hcat(samples.adc.t, samples.adc.A)
    seq_dict = Dict(
        "Gx" => Gx,
        "Gy" => Gy,
        "Gz" => Gz,
        "RF_AM" => RF_AM,
        "RF_FM" => RF_FM,
        "ADCS" => ADCS
    )

    # Write to matlab file
    matwrite(joinpath(matfolder, matfilename), Dict("sequence" => seq_dict))
end

function export_2_mat_kspace(seq, matfolder; matfilename="seq_kspace.mat")
    sampling_rule = MaxStepSizeRule(1, 5e-5)
    kspace, kspace_adc = get_kspace(seq; sampling_rule)
    matwrite(joinpath(matfolder, matfilename), Dict("kspace" => kspace, "kspace_adc" => kspace_adc))
end

function export_2_mat_moments(seq, matfolder; matfilename="seq_moments.mat")
    sampling_rule = MaxStepSizeRule(1, 5e-5)
    seqd = discretize(seq; sampling_rule)
    t = seqd.t[1:end-1]
    k0, _ =  KomaMRIBase.get_kspace(seqd)
    k1, _ =  KomaMRIBase.get_M1(seqd)
    k2, _ =  KomaMRIBase.get_M2(seqd)
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
    sys_dict = Dict("B0" => sys.limits.B0,
                "B1" => sys.limits.B1,
                "Gmax" => sys.limits.Gmax,
                "Smax" => sys.limits.Smax,
                "ADC_dt" => sys.limits.ADC_Δt,
                "DUR_dt" => sys.limits.DUR_Δt,
                "GR_dt" => sys.limits.GR_Δt,
                "RF_dt" => sys.limits.RF_Δt,
                "RF_ring_down_time" => sys.limits.RF_ring_down_time,
                "RF_dead_time" => sys.limits.RF_dead_time,
                "ADC_dead_time" => sys.limits.ADC_dead_time)
    matwrite(joinpath(matfolder, matfilename), Dict("scanner" => sys_dict))
end

function export_2_mat_raw(raw_ismrmrd, matfolder; matfilename="raw.mat");
    if haskey(raw_ismrmrd.params, "userParameters")
        dict_for_mat = Dict()
        dict_user_params = raw_ismrmrd.params["userParameters"]
        for (key, value) in dict_user_params
            ascii_key = replace(key, "Δ" => "d") # MATLAB does not support non-ascii characters in variable names
            dict_for_mat[ascii_key] = value
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
    files = String[]
    if type=="all"
        export_2_mat_sequence(seq, matfolder)
        export_2_mat_kspace(seq, matfolder)
        export_2_mat_moments(seq, matfolder)
        export_2_mat_phantom(phantom, matfolder)
        export_2_mat_scanner(sys, matfolder)
        export_2_mat_raw(raw_ismrmrd, matfolder)
        export_2_mat_image(image, rec_params, matfolder)
    elseif type=="sequence"
		files = [head*"_sequence.mat", head*"_kspace.mat", head*"_moments.mat"]
		export_2_mat_sequence(seq, matfolder; matfilename=files[1])
        export_2_mat_kspace(seq, matfolder; matfilename=files[2])
        export_2_mat_moments(seq, matfolder; matfilename=files[3])
	elseif type=="phantom"
		files = [head*"_phantom.mat"]
		export_2_mat_phantom(phantom, matfolder; matfilename=only(files))
    elseif type=="scanner"
		files = [head*"_scanner.mat"]
		export_2_mat_scanner(sys, matfolder; matfilename=only(files))
    elseif type=="raw"
		files = haskey(raw_ismrmrd.params, "userParameters") ? [head*"_raw.mat", "sim_params.mat"] : String[]
		isempty(files) || export_2_mat_raw(raw_ismrmrd, matfolder; matfilename=first(files))
    elseif type=="image"
		files = [head*"_image.mat"]
		haskey(rec_params, :reconSize) && push!(files, "rec_params.mat")
		export_2_mat_image(image, rec_params, matfolder; matfilename=first(files))
	end

    isempty(files) && return "<ul><li><b>Path:</b> $matfolder</li></ul>"
    label = length(files) == 1 ? "Name" : "Names"
    return "<ul><li><b>$label:</b> $(join(files, ", "))</li><li><b>Path:</b> $matfolder</li></ul>"
end
