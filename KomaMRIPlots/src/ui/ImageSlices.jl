function image_slice_plot(image; spacing=(1,1,1), zmin=minimum(image), zmax=maximum(image), kwargs...)
    spacing = Float64.(collect(spacing) ./ first(spacing))
    shape = size(image)
    positions = map(n -> collect(0:n-1), shape)
    center = collect(shape .÷ 2)
    planes = map(PHANTOM_SLICE_PLANES) do plane
        (; plane..., name=join((('x','y','z')[axis] for axis in plane.dimensions), "-"))
    end
    plot = plot_image(zeros(2,2); zmin,zmax,kwargs...)
    layout = plot.layout
    template = only(plot.data)
    slice_layout!(layout, planes, positions, center, map(n -> [-0.5,n-0.5],shape),spacing) do axis,value
        "$(('x','y','z')[axis]) = $value"
    end
    layout[:margin] = attr(; l=20,r=20,t=30,b=38,autoexpand=false)
    layout.fields[:meta] = Dict{String,Any}("image_slices" => (; planes,center,bounds=map(n -> [0,n-1],shape),
        shape,spacing,values=[vec(image)],zlabels=positions[3]))
    layout[:uirevision] = "image-slices-$(join(shape,'-'))"
    plot.data = map(enumerate(planes)) do (panel,plane)
        trace = deepcopy(template)
        horizontal,vertical = plane.dimensions
        trace[:x],trace[:y] = positions[horizontal],positions[vertical]
        trace[:xaxis],trace[:yaxis] = plane.refs
        trace[:z] = [[image[ntuple(axis -> axis==horizontal ? x+1 : axis==vertical ? y+1 : center[axis]+1,3)...]
            for x in positions[horizontal]] for y in positions[vertical]]
        trace[:customdata] = [[collect(ntuple(axis -> axis==horizontal ? x : axis==vertical ? y : center[axis],3))
            for x in positions[horizontal]] for y in positions[vertical]]
        trace[:showscale] = panel==3
        trace[:colorbar] = attr(; x=0.76,xanchor="center",y=0.2325,len=0.465)
        trace[:hovertemplate] = "x: %{customdata[0]}<br>y: %{customdata[1]}<br>z: %{customdata[2]}<br>Value: %{z:.4g}<extra></extra>"
        trace
    end
    return plot
end

plot_image(image::AbstractArray{<:Real,3}; kwargs...) = image_slice_plot(image; kwargs...)
