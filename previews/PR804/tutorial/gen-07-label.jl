using KomaMRI, PlotlyJS, Suppressor #hide
sys = Scanner(); #hide
Nx = 64
Ny = 64
Nslices = 2
FOV = 23e-2 #hide
Tadc = sys.ADC_Δt * (Nx - 1) #hide
TR = 0.45 #hide
slice_delay = 2.0 #hide
TE = 5e-3 #hide
Gx = 1 / (γ * sys.ADC_Δt * FOV) #hide
Gz = 20e-3 #hide
Trf = 2e-3 #hide
slice_positions = range(-4e-3, 4e-3; length=Nslices); #hide
obj = brain_phantom3D(); #hide
obj.Δw .= 0; #hide
phase_encode(lin) = Grad((lin - (Ny - 1) / 2) * Gx * sys.ADC_Δt / Tadc, Tadc) #hide
rf = [PulseDesigner.RF_sinc(sys.B1, Trf, sys; G=[0, 0, Gz], Δf=γ * Gz * z) for z in slice_positions]; #hide

function cartesian_label_sequence()
    seq = Sequence()
    adc = ADC(Nx, Tadc)
    readout = Grad(Gx, Tadc)
    prewinder = Grad(-Gx / 2, Tadc)
    @addblocks for slc in 0:(Nslices - 1), lin in 0:(Ny - 1)
        pe = phase_encode(lin)
        seq += rf[slc + 1]
        seq += (x=prewinder, y=pe)
        seq += Delay(TE) #hide
        seq += (x=readout, adc, LabelSet(slc, "SLC"), LabelSet(lin, "LIN"))
        seq += (y=-pe) #hide
        seq += Delay(TR)
        if lin == Ny - 1 #hide
            seq += Delay(slice_delay) #hide
        end #hide
    end
    return seq
end

seq = cartesian_label_sequence();

raw = @suppress simulate(obj, seq, sys)

for p in [1, Ny, Ny + 1, Nslices * Ny]
    slc = raw.profiles[p].head.idx.slice |> Int
    lin = raw.profiles[p].head.idx.kspace_encode_step_1 |> Int
    println("Profile $p: SLC=$slc, LIN=$lin")
end

raw.params["trajectory"] = "cartesian";
raw.params["encodedSize"] = [Nx, Ny, 1];
raw.params["reconSize"] = [Nx, Ny, 1];

acqData = AcquisitionData(raw; estimateProfileCenter=true);
rec = reconstruction(acqData, Dict{Symbol,Any}());
image = abs.(rec);
image_display = sqrt.(image ./ maximum(image)); #hide

p1 = plot_image(image_display[:, :, 1]; height=360, title="Slice 1", zmin=0, zmax=1)
p2 = plot_image(image_display[:, :, 2]; height=360, title="Slice 2", zmin=0, zmax=1)
foreach(p -> p.plot.data[1].showscale = false, (p1, p2)); #hide
display([p1 p2]);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
