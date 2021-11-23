path = @__DIR__

# loadbutton = filepicker()
# columnbuttons = Observable{Any}(dom"div"())
# data = Observable{Any}(Phantom)
# map!(f->begin
#         global recon = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"recon")
#         print("Loaded! ... $f\n")
#     end
#     , data, loadbutton)

# Ploting recon
global recon = ifftc(kdata)

hh, ww = 420,550
# l = PlotlyJS.Layout(;title="Reconstruction", yaxis_title="y", yaxis=attr(scaleanchor="x"),
#     xaxis_title="x",height=hh,width=ww,
#     modebar=attr(orientation="v"),scene=attr(aspectratio=attr(x=1,y=1,z=1)))
l = PlotlyJS.Layout(;title="Reconstruction",yaxis_title="y",
    xaxis_title="x",height=hh,width=ww,
    yaxis=attr(scaleanchor="x"),
    modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
p = PlotlyJS.plot(PlotlyJS.heatmap(z=abs.(recon),showscale=false,colorscale="Greys",transpose = false),l)
plt = Observable{Any}(p)
# PlotlyJS.savefig(p, path*"/assets/phantom.png", width=320, height=300)

recon = dom"div"(plt)
content!(w, "div#content", recon)
