using KomaMRI, PlotlyJS

seq_file = "./test/ourSequence/gre3d.seq" #"./examples/1.sequences/ge.seq" 
GRE = read_seq(seq_file)
GRE.DEF = Dict("Name"=>"GRE3D", "Nx" => 64, "Ny" => 64, "Nz" => 32)
obj = brain_phantom3D()
sys = Scanner()
sys.B0 = 3
sims = Dict{String, Any}("Δt" => 0.001, "Δt_rf" => 6.4e-6) 
#savefig(plot_seqd(GRE, simParams=sims),"./test/ourSequence/DiscreteSeqMid.png")
rawfile = "./test/ourSequence/Koma_signal.mrd"

#savefig(plot_seq(GRE),"./test/ourSequence/Sequence.png")
#savefig(plot_seq(GRE,range=[2200,2280]),"./test/ourSequence/SequenceMed.png")
#savefig(plot_seq(GRE,range=[0,2000]),"./test/ourSequence/SequenceEarly.png")

raw_ismrmrd = simulate(obj,GRE,sys,simParams=sims)
show(IOContext(stdout, :limit => true), "text/plain", raw_ismrmrd)
p3 = plot_signal(raw_ismrmrd; slider=false, height=300)
savefig(p3, "./test/ourSequence/Signal.png")
fout = ISMRMRDFile(rawfile)
#KomaMRICore.save(fout, raw_ismrmrd)
acq = AcquisitionData(raw_ismrmrd)
Nx, Ny, Nz = raw_ismrmrd.params["reconSize"][1:3]
reconParams = Dict{Symbol,Any}(:reco=>"standard", :reconSize=>(Ny,Nx,Nz))
image = reconstruction(acq, reconParams)
for n in 1:size(image,3)
    local Title = "Slice " * string(n)
    local Save = "./test/ourSequence/slices/StrangSplitting/20 FOV 64x64x32/Slice " * string(n) * ".png"
    savefig(plot_image(abs.(image[:, :, n]); height=360, title=Title),Save)
end
