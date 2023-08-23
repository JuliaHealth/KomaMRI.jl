using KomaMRI,DelimitedFiles
Fields = readdlm("./test/ourSequence/mpf_2023062211393219.dat",' ')
Gx = Float16.(Fields[:,1] .* 10^(-3))
Gy = Float16.(Fields[:,2] .* 10^(-3))
Gz = Float16.(Fields[:,3] .* 10^(-3))
Rf = Float32.(Fields[:,4] .* 10^(-6))   
N = length(Gx)
dt = 6.4 * 10^(-6)
println(typeof(Gx))
GRE = Sequence()
println(N)
replace!(Rf, NaN => 0)
replace!(Gx, NaN => 0)
replace!(Gy, NaN => 0)
replace!(Gz, NaN => 0)
i = 46754
#
for n in 1:3
    global i
    global GRE
    local Grads = [Grad(Gx[i+87],308*dt,87*dt) Grad(0,0) Grad(Gx[i+837],0,55*dt) Grad(Gx[i+941],80*dt,50*dt) Grad(Gx[i+941+ 80 + 50*2],80*dt,50*dt) Grad(Gx[i+1469],476*dt,217*dt,88*dt,0);
             Grad(0,0) Grad(0,0) Grad(0,0) Grad(0,0) Grad(0,0) Grad(0,0);
             Grad(0,0) Grad(Gz[i+484],284*dt,7*dt) Grad(Gz[i+837],0,55*dt)  Grad(0,0) Grad(0,0) Grad(0,0);]
    i+= 31250	
    GRE += Sequence(Grads)
    if (n != 3)
        GRE += Sequence([Grad(0,31250*dt);;],[RF(0,0);;],[ADC(100,31250*dt)])
    end
end		
i = 114368
prevDif = 3249
for n in 1:3
    global i
    global GRE
    local Rfs = [RF(Rf[i],29*dt,0,prevDif*dt) ;;]
    global prevDif = 595
    GRE += Sequence([Grad(0,0);;],Rfs)
    i+= 596 + 29
end	
GRE += Sequence([Grad(0,0);;],[RF(0,0);;],[ADC(100,77937*dt)])
i = 194211
Grads = [Grad(Gx[i+ 968],0,235*dt,733*dt) Grad(Gx[i + 2209],0,235*dt,118*dt) Grad(Gx[i + 4249],0,235*dt,118*dt);
    Grad(0,0) Grad(0,0) Grad(0,0); 
    Grad(Gz[i+233],1623*dt,233*dt,0,0) Grad(Gz[i + 233*2 + 1623],2041*dt,0,0,0) Grad(Gz[i + 233*2 + 1623 + 2041], 1722*dt,0,231*dt,0)]
Rfs = [RF(Rf[i + 617:i + 730],113*dt,0,617*dt) RF(Rf[i + 1858:i + 1971],113*dt) RF(Rf[i + 3899:i + 4012],113*dt);]
#Adc = [ADC(100, 617*dt), ADC(100,1928*dt,113*dt), ADC(100,1953*dt,113*dt)]
GRE += Sequence(Grads,Rfs) 
i = 277984
GRE += Sequence([Grad(0,0,0,0,77919*dt) ;;],[RF(0,0);;],[ADC(100,77919*dt)])
Grads = [Grad(Gx[i+ 968],0,235*dt,733*dt) Grad(Gx[i + 2209],0,235*dt,118*dt) Grad(Gx[i + 4249],0,235*dt,118*dt) Grad(Gx[i + 6808],416*dt,327*dt,953*dt);
    Grad(0,0) Grad(0,0) Grad(0,0) Grad(Gy[i + 6808],416*dt,327*dt,953*dt); 
    Grad(Gz[i+233],1623*dt,233*dt,0,0) Grad(Gz[i + 233*2 + 1623],2041*dt,0,0,0) Grad(Gz[i + 233*2 + 1623 + 2041], 1722*dt,0,231*dt,0) Grad(Gz[i + 6808],416*dt,327*dt,953*dt)]
Rfs = [RF(Rf[i + 617:i + 730],113*dt,0,617*dt) RF(Rf[i + 1858:i + 1971],113*dt) RF(Rf[i + 3899:i + 4012],113*dt) RF(0,0);]
GRE += Sequence(Grads,Rfs)

# !!! Important stuff at line 5510631

GRE += Sequence([Grad(0,0,0,0,49004*dt);;],[RF(0,0);;],[ADC(100,49004*dt)])
i = 334539
while (i + 2*5518) < 5510100 # cannot run the entire thing  
    global i
    global GRE
    #local Grads = [Grad(0,0) Grad(Gx[i+596+164],403*dt,(164)*dt) Grad(0,0) Grad(Gx[i+2329],0,45*dt,235*dt,0);
    #        Grad(0,0) Grad(0,0) Grad(0,0) Grad(Gy[i+2347],160*dt,60*dt);
    #        Grad(Gz[i+8],581*dt, 8*dt) Grad(Gz[i+760],406*dt,164*dt) Grad(Gz[i+1561],723*dt,231*dt,0,0) Grad(Gz[i+2347],0,60*dt,220*dt,0);] 
    local Grads = [Grad(0,0) Grad(Gx[i+596+164],403*dt,(164)*dt) Grad(0,0) Grad(Gx[i+2329],0,45*dt,235*dt,0);
            Grad(0,0) Grad(0,0) Grad(0,0) Grad(Gy[i+2347],160*dt,60*dt);
            Grad(Gz[i+8],581*dt, 8*dt) Grad(Gz[i+760],406*dt,164*dt) Grad(Gz[i+1561],723*dt,231*dt,0,0) Grad(Gz[i+2347],0,60*dt,220*dt,0);]       
    local Rfs = [RF(Rf[i + 93:i + 587],495*dt, 0 ,93*dt) RF(0,0) RF(Rf[i+1565:i+2280],716*dt, 0,234*dt) RF(0,0);]		
    local Adc = [ADC(5,324*dt,50*dt), ADC(5, 335*dt, 125*dt), ADC(5, 335*dt, 125*dt), ADC(5, 335*dt, 125*dt), ADC(5, 335*dt,125*dt)]
    GRE += Sequence(Grads,Rfs)
    xGrad = [Grad(Gx[i+2578],324*dt +50*dt+115*dt) Grad(Gx[i+3152],335*dt+125*dt*2) Grad(Gx[i+3152+335+125*2],335*dt + 125*dt*2) Grad(Gx[i+3152+335*2+125*4],335*dt + 125*dt*2) Grad(Gx[i+3152+335*3+125*6],335*dt + 125*dt + 256*dt,0);]
    GRE += Sequence(xGrad,[RF(0,0) RF(0,0) RF(0,0) RF(0,0) RF(0,0);])
    #GRE += Sequence([Grad(0,2961*dt) ;;])
    i+= 5519
    local Grads = [Grad(0,0) Grad(Gx[i+596+164],403*dt,(164)*dt) Grad(0,0) Grad(Gx[i+2329],0,45*dt,235*dt,0);
            Grad(0,0) Grad(Gy[i+596+164],407*dt,164*dt) Grad(0,0) Grad(Gy[i+2347],160*dt,60*dt);
            Grad(Gz[i+8],581*dt, 8*dt) Grad(Gz[i+760],406*dt,164*dt) Grad(Gz[i+1561],723*dt,231*dt,0,0) Grad(Gz[i+2347],0,60*dt,220*dt,0);]
    local Rfs = [RF(Rf[i + 93:i + 587],495*dt, 0 ,93*dt) RF(0,0) RF(Rf[i+1565:i+2280],716*dt, 0,234*dt) RF(0,0);]	
    GRE += Sequence(Grads,Rfs)
    GRE += Sequence(xGrad,[RF(0,0) RF(0,0) RF(0,0) RF(0,0) RF(0,0);],Adc)
    i+= 5519
end	
#i = 5510171
GRE.DEF = Dict("Name"=>"GRE3D", "FOV" => [230, 160.178574, 96])
obj = brain_phantom3D()
obj.Δw .= 0
sys = Scanner()
sims = Dict{String, Any}("Δt" => 0.001, "Δt_rf" => 6.4e-6) 
println("Hello")
#savefig(plot_seqd(GRE, simParams=sims),"./test/ourSequence/DiscreteSeqMid.png")
rawfile = "./Koma_signal.mrd"
#savefig(plot_seq(GRE),"./test/ourSequence/Sequence.png")
#savefig(plot_seq(GRE,range=[2200,2280]),"./test/ourSequence/SequenceMed.png")
#savefig(plot_seq(GRE,range=[0,2000]),"./test/ourSequence/SequenceEarly.png")
raw_ismrmrd = simulate(obj,GRE,sys,simParams=sims)
show(IOContext(stdout, :limit => true), "text/plain", raw_ismrmrd)
p3 = plot_signal(raw_ismrmrd; slider=false, height=300)
savefig(p3, "./test/ourSequence/Signal.png")
fout = ISMRMRDFile(rawfile)
KomaMRICore.save(fout, raw_ismrmrd)
