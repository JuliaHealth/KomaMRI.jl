# This is a test of the MRIReco package. In the future I plan to change the reconstruction step
# to use this. I will need to change the raw data output to a ISMRMRDFile.

using MRIReco
 
f = ISMRMRDFile("/home/ccp/Downloads/knee_3dFSE_slice170.h5")
acqData = AcquisitionData(f);

recoParams = Dict{Symbol, Any}()
recoParams[:reco] = "direct"
recoParams[:reconSize] = (320,320) # this size is also contained in acqData.encodingSize
img = reconstruction(acqData, recoParams)
# #samplingDensity.jl, v \ u uses too much memory

# ##
img = sqrt.(sum(img.^2,dims=5))
using Plots
heatmap(reverse( abs.(img[:,:,1,1,1]),dims=1), c=:viridis)

# multiCoil
# sampling
redFac= 4.0
acqDataSub = sample_kspace(acqData,redFac,"poisson",calsize=30,profiles=false);
# show sampling pattern
msk = zeros(acqDataSub.encodingSize[1],acqDataSub.encodingSize[2])
msk[acqDataSub.subsampleIndices[1]] .= 1
# #sampling 
# heatmap(reverse( msk,dims=1), c=:grays)

#ESPIRIT
# coil sensitivities
smaps = espirit(acqData,(6,6),30,eigThresh_1=0.035,eigThresh_2=0.98)
# SENSE reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320,320) #previously it says :shape
params[:senseMaps] = smaps
params[:solver] = "cgnr"
params[:regularization] = "L2"
params[:λ] = 1e-4
params[:iterations] = 1
params[:normalizeReg] = true

img_reco = reconstruction(acqDataSub, params)

## show image
heatmap(reverse( abs.(img_reco[:,:,1,1,1]),dims=1), c=:viridis)

##
using MRIsim

imcoils = []
Nc = numChannels(data) 
for coil = 1:Nc
    raw = reshape(data.kdata[1][:,coil], sz[1:2]...)
    img = MRIsim.ifftc(raw)
    append!(imcoils, [img])
end
imgss = sqrt.( sum([abs.(imcoils[i]).^2 for i=1:Nc]) ) #sum of squares

##

l = PlotlyJS.Layout(;height=500, width=500, colorbar=false)
p = PlotlyJS.heatmap(z = abs.(imgss), colorscale="Greys")
PlotlyJS.plot(p,l)
