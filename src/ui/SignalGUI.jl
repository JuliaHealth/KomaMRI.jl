# loadbutton = filepicker()
# columnbuttons = Observable{Any}(dom"div"())
# data = Observable{Any}(Array{Sequence})
# map!(f->begin
#         print("Loading... $f\n")
#         JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
#     end
#     , data, loadbutton)
Tmin = minimum(t_interp)*1e3
Tmax = maximum(t_interp)*1e3
Tacq = Tmax - Tmin
hh, ww = 450, 150
l = PlotlyJS.Layout(;title="Acquired signal", yaxis_title="Signal (a.u.)",
    xaxis_title="time [ms]",
    xaxis=attr(range=[Tmin+Tacq/2-.025*Tacq,Tmin+Tacq/2+.025*Tacq],
			    rangeslider=attr(visible=true),
                gridcolor="gray"),
    yaxis=attr(gridcolor="gray"),
    plot_bgcolor="black",
    modebar=attr(orientation="v"),legend=false,height=400)
absS = PlotlyJS.scatter(x=t_interp*1e3,y=abs.(signal), name="|S(t)|")
reS =  PlotlyJS.scatter(x=t_interp*1e3,y=real.(signal),name="Re{S(t)}")
imS =  PlotlyJS.scatter(x=t_interp*1e3,y=imag.(signal),name="Im{S(t)}")
p = PlotlyJS.plot([absS,reS,imS],l)
plt = Observable{Any}(p)
# function makebuttons(df)
#     namesdf = 1:length(df)
#     buttons = button.(string.())
#     for (btn, name) in zip(buttons, namesdf)
#         map!(t -> begin
#             p = PlotlyJS.plot(df[name][:])
#             PlotlyJS.relayout!(p,height=hh)
#             p
#             end
#             , plt, btn)
#     end
#     dom"div"(hbox(buttons))
# end
# map!(makebuttons, columnbuttons, data)
ui = dom"div"(plt)
content!(w, "div#content", ui)
