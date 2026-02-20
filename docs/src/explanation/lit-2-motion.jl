# # Motion

using KomaMRI #hide
using PlotlyJS, Random #hide
obj = brain_phantom2D(); #hide

# Koma can easily simulate the effects of motion during acquisitions. 
# As introduced in the [previous section](1-phantom.md), the motion-related information
# of the phantom is stored in the `motion` field of its structure.

# Koma's motion model has been designed to accomodate a variety of real-world scenarios, including:

# - [Patient motion inside a scanner](../tutorial/05-SimpleMotion.md), which may involve simultaneous or sequential translations and rotations of body parts during the acquisition.
# - Myocardial motion, including simulataneous contraction, rotation, torsion, and translation motion within the cardiac cycle.
# - Pseudo-periodic heart patterns, caused by variations in heart rate or arrhythmias that prevent the heart's motion from being perfectly periodic.
# - Flow through blood vessels, where the spin trajectories or fluid fields may have been obtained from Computational Fluid Dynamics (CFD) simulations.
# - [Diffusion](../tutorial/06-DiffusionMotion.md), which can be modeled, among many other ways, as microscopic Brownian spin trajectories.

# ... And, ultimately, any type of motion you can think of, no matter how complex!

# To handle these scenarios, Koma represents motion as a collection of elementary movements that can be independently 
# configured and combined. This approach allows for the definition of any complex motion pattern, with the ability to 
# specify overlapping time intervals and even model bidirectional motions along predefined trajectories.

# ## Understanding the `motion` field and its possible values
# The `motion` field within the `Phantom` struct can take different values depending on whether the phantom is static or dynamic. 
# For static phantoms, the field is set to `NoMotion`. For dynamic phantoms, the field can be either a `Motion` or a `MotionList` struct.
# A `Motion` represents a single movement, characterized by an `action`, a `time` curve, and a range of affected `spins`. 
# Regarding the `MotionList` struct, it is simply a collection of `Motion` instances, which is useful for defining motion compositions.

# ```julia
# struct Phantom{T<:Real}
#     (...)
#     #Motion
#     motion::Union{NoMotion, Motion{T}, MotionList{T}} = NoMotion()
# end
# ```

# ### `NoMotion` struct
# `NoMotion` is the default type for static phantoms. Since its structure has no fields, making a phantom static is as simple as:
obj.motion = NoMotion();

# ### `Motion` struct
# The `Motion` struct contains information about a basic motion, understood as the combination of an `action`, a `time` curve and a `spins` span.
# This three fields will be described in detail later. Here is an example of how to assign a motion to a phantom in this case:
obj.motion = Motion(Translate(0.0, 0.1, 0.2), TimeRange(0.0, 1.0), AllSpins());

# !!! note
#     There are `Motion` constructors that simplify its definition and have the same name as the actions, but written in lowercase. For example:
# 
#     ```julia
#     obj.motion = translate(0.0, 0.1, 0.2, TimeRange(0.0, 1.0), AllSpins())
#     ```
#
#     This is equivalent to writing:
#     ```julia
#     obj.motion = Motion(Translate(0.0, 0.1, 0.2), TimeRange(0.0, 1.0), AllSpins())
#     ```

# ### `MotionList` struct
# The `MotionList` struct contains a single field called `motions`, which is a vector of `Motion` instances.
# This design makes it possible to define both sequential and simultaneous concatenations of motions over time.
# An example of how this would be used is:
obj.motion = MotionList(
    Motion(Translate(0.0, 0.1, 0.2), TimeRange(0.0, 1.0), AllSpins()),
    Motion(Rotate(0.0, 0.0, 45.0, (0.0, 0.0, 0.0)), Periodic(1.0, 0.5), SpinRange(1:1000))
);

# ## The `Motion` structure and its fields
# The `Motion` struct is the basic building block for defining motion in Koma. 
# As we mentioned earlier, it has three main fields: `action`, `time`, and `spins`. 
# Together, these fields define what the motion is, when it happens, and which spins are involved:
# ```julia
# struct Motion{T<:Real}
#     action::AbstractAction{T}
#     time  ::TimeCurve{T}
#     spins ::AbstractSpinSpan
# end
# ```

# ### The `action` field
# Let's start with the `action` field, which defines the type and magnitude (i.e., the final state) of the motion.
# Currently, Koma supports five actions: [`Translate`](@ref), [`Rotate`](@ref), [`HeartBeat`](@ref), [`Path`](@ref), and [`FlowPath`](@ref).
# The first three fall under the category of _SimpleActions_, while the last two belong to the _ArbitraryActions_.
# _SimpleActions_ are defined by parameters that are easy to understand and use, such as translation distance, rotation angles, or contraction rates.
# _ArbitraryActions_, on the other hand, are more complex and can be defined by a set of spin trajectories.

# ### The `time` field
# The `time` field defines how the motion behaves over time and must be an instance of the `TimeCurve` struct, which
# works similarly to animation curves in video editing, 3D design, or video games. Essentially, it allows you to adjust the "timing" of 
# the motion without affecting its magnitude or other characteristics.

# Given an initial and final state (see the `action` field), time curves allow you to define how the transition between 
# those states should occur. The `TimeCurve` structure lets you define an animation curve by specifying the coordinates
# of its points, along with two additional parameters that control its periodicity and pseudo-periodicity:

# ```julia
# struct TimeCurve{T<:Real}
#     t::AbstractVector{T}
#     t_unit::AbstractVector{T}
#     periodic::Bool                 
#     periods::Union{T,AbstractVector{T}}
# end
# ```

# This enables you to create any type of curve, and thus, any kind of motion pattern over time.

# A full description of this structure, including examples and constructors, can be found in the [`TimeCurve` API reference](../reference/2-koma-base.md#KomaMRIBase.TimeCurve).

# ### The `spins` field
# Finally, the `spins` field must be an instance of the [`AbstractSpinSpan`](../reference/2-koma-base.md#AbstractSpinSpan-types) type.
# It defines which spins in the phantom are affected by the motion, and which of them remain static.
# This allows you to define motions that only affect a subset of spins, while keeping others unaffected.

# ## See it in action
# Now that we have a basic understanding of the `motion` field and its components, let's see some usage examples.
# In all cases, we start with the same phantom: a hollow cube with 1 mm side length and 20 µm spin spacing,
# centered at the origin and aligned with the coordinate axes. To make the motion easier to visualize, 
# each face of the cube is given a different T1 value:

#jl # Phantom construction
L = 1e-3 #hide
Δx = 20e-6 #hide
x = y = z = -L/2:Δx:L/2 #hide
xx = reshape(x, (length(x),1,1)) #hide
yy = reshape(y, (1,length(y),1)) #hide 
zz = reshape(z, (1,1,length(z))) #hide 
x = 1*xx .+ 0*yy .+ 0*zz #hide
y = 0*xx .+ 1*yy .+ 0*zz #hide
z = 0*xx .+ 0*yy .+ 1*zz #hide
◼(L) =((abs.(x) .<= L/2) .& (abs.(y) .<= L/2) .& (abs.(z) .<= L/2)) #hide
cube = ◼(L) - ◼(L - Δx) # Hollow cube #hide 
ρ = 1.0*cube #proton density #hide
T1 = copy(ρ) #hide
T1s = [100, 500, 1000, 2500, 2000, 1500] .* 1e-3 #hide
idx_T1 = 1 #hide
ϵ = 1e-5 #hide
for (i, x) in enumerate([x,y,z]) #hide
    for (j, L) in enumerate([-L/2, L/2]) #hide
        T1[(L - ϵ) .<= x .<= (L + ϵ)] .= T1s[idx_T1] #hide
        global idx_T1 += 1 #hide
    end #hide
end #hide
obj = Phantom( x=x[ρ .!= 0], y=y[ρ .!= 0], z=z[ρ .!= 0], T1 = T1[ρ .!= 0] ); #hide
p = plot_phantom_map(obj, :T1; height=440); #hide
#jl display(p);

# ### Translation motion
# In this first example, we've added a translational motion of -0.5, 0.6, and 0.7 mm along the
# three spatial directions. The motion lasts for 1 second and affects the entire phantom:

obj.motion = translate(-5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), AllSpins());

# Let’s plot this phantom and see how it moves. The `time_samples` argument specifies the number of time samples to be plotted.
# You can use the bottom slider to scroll through time and check its exact position at each moment: 

p1 = plot_phantom_map(obj, :T1; time_samples=11, height=440)
#jl display(p1)

# ### Rotation motion
# In this case, we add a rotational motion to the phantom: 90º around the y-axis and 75º around
# the z-axis. Like before, the motion lasts for 1 second and affects all spins in the phantom:

obj.motion = rotate(0.0, 90.0, 75.0, TimeRange(0.0, 1.0), AllSpins());

p2 = plot_phantom_map(obj, :T1; time_samples=11, height=440) #hide
#jl display(p2);

# ### Adding motion to a phantom subset
# Sometimes, you may want to assign motion to just a part of the phantom instead of the whole thing. 
# This can be done using the [`SpinRange`](@ref) structure, where you specify the indices of the spins that
# should be affected. In this example, we apply a translational motion to the upper half of the phantom:

obj.motion = translate(-5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), SpinRange(7500:15002));

p3 = plot_phantom_map(obj, :T1; time_samples=11, height=440); #hide
#jl display(p3);

# ### Motion combination
# You can freely add multiple motions to a phantom, each with its own type, time span, 
# and affected spin range. These motions can overlap in time (affecting the phantom simultaneously) 
# or happen one after another. Both cases are fully supported, so you're free to combine different 
# effects across various parts of the phantom and time intervals, creating as complex a motion pattern as you need.

# This example shows two brain phantoms undergoing the same translational and rotational motions, but with
# different time spans. In the top phantom, the translation takes place from 0 to 0.5 seconds, followed by the 
# rotation from 0.5 to 1 second. In the bottom phantom, both motions happen over the same time span, from 0 to 1 second:

obj1 = brain_phantom2D() #hide
obj2 = copy(obj1) #hide
obj1.x .-= 20e-2; obj2.x .-= 20e-2 #hide
obj1.y .+= 12e-2; obj2.y .-= 12e-2 #hide
obj1.motion = MotionList(
    translate(40e-2, 0.0, 0.0, TimeRange(0.0, 0.5),AllSpins()),
    rotate(0.0, 0.0, 90.0, TimeRange(0.5, 1.0),AllSpins()),
)

obj2.motion = MotionList(
    translate(40e-2, 0.0, 0.0, TimeRange(0.0, 1.0),AllSpins()),
    rotate(0.0, 0.0, 90.0, TimeRange(0.0, 1.0),AllSpins()),
)

obj = obj1 + obj2
p4 = plot_phantom_map(obj, :T1; time_samples=11, view_2d=true, height=440) #hide
#jl display(p4);

# ### Realistic head motion
# As a more realistic final example, let's try to replicate the head motion made by a patient inside the scanner.
# This motion consists of a series of translations and rotations, with the rotation center being the neck:
Random.seed!(1234) #hide  
obj = brain_phantom2D()

Nintervals = 10  
interval_dur = 0.1 
tra_x = rand(-5:5, Nintervals) .* 1e-3 # Translation in x [m]  
tra_y = rand(-5:5, Nintervals) .* 1e-3 # Translation in y [m]  
rot_z = rand(-5:5, Nintervals) .* 1e0  # Rotation in z    [º]  
rot_center = (0.0, -3.0, 0.0)  .* 1e-2 # Rotation around the neck  

motion_list = Motion[] 
for i in 1:Nintervals  
    t_interval = TimeRange(interval_dur * (i-1), interval_dur * i)  
    tra = translate(tra_x[i], tra_y[i], 0.0, t_interval) 
    rot = rotate(0.0, 0.0, rot_z[i], t_interval; center=rot_center)  
    push!(motion_list, [tra, rot]...)  
end  

obj.motion = MotionList(motion_list...);  

p5 = plot_phantom_map(obj, :T1; time_samples=21, view_2d=true, height=440) #hide
#jl display(p5);

# The motion signals for this phantom are shown in the plot below, where you can see the translations in x and y, and the rotation around z over time.
p6 = plot( #hide
    (0:interval_dur:interval_dur*length(tra_x)) .* 1e3, #hide
    [cumsum([0, tra_x...]) * 1e3 cumsum([0, tra_y...]) * 1e3 cumsum([0, rot_z...])], #hide
    Layout( #hide
        title = "Head motion profile", #hide
        xaxis_title = "time (ms)", #hide
        yaxis_title = "Position" #hide
    )) #hide
restyle!(p6,1:3, name=["X-Trans (mm)", "Y-Trans (mm)", "Z-Rot (º)"]) #hide
#md p6 #hide
#jl display(p6);

# A simulation and motion-corrected reconstruction based on a similar, slightly simplified head motion is available [here](../tutorial/05-SimpleMotion.md).
