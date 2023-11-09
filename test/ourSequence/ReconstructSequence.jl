using MRIReco
using KomaMRI

f = ISMRMRDFile("./test/ourSequence/Koma_signal.mrd")

raw = RawAcquisitionData(f)
acq = AcquisitionData(raw)

objDict = Dict{String,Any}()
objDict["raw"] = raw
objDict["acq"] = acq

println(acq.traj[1])
acq.traj[1].numProfiles = 2048

temp_prof = acq.kdata[1]

acq.encodingSize = (64,64,32)

acq.traj[1].nodes[1,:] *= 5.0
acq.traj[1].nodes[2,:] *= 4.98

Nx = 64
Ny = 64
Nz = 32

reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx,Ny,Nz))
image = reconstruction(acq, reconParams)

for n in 1:size(image,3)
    local Title = "Slice " * string(n)
    local Save = "./test/ourSequence/slices/StrangSplitting/20 FOV 64x64x32/Slice " * string(n) * ".png"
    savefig(plot_image(abs.(image[:, :, n]); height=360, title=Title),Save)
end


