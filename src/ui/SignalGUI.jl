#plt = Observable{Any}(plot_signal(raw_ismrmrd; darkmode))
#ui = dom"div"(plt)
#map!(p->plot_signal(p; darkmode), plt, sig_obs)
#content!(w, "div#content", ui)

plt = Observable{Any}(plot_signal(raw_ismrmrd; darkmode))
btn = button("Export .mat")
#ui = dom"div"(vbox(dom"div"(hbox(btn)), plt))
ui = dom"div"(plt)
map!(p->plot_signal(p; darkmode), plt, sig_obs)
content!(w, "div#content", ui)

function export2mat()

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

    matwrite("raw.mat", Dict("raw" => raw_dict))

end

on(n -> export2mat(), btn)
