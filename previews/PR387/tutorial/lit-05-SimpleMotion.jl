# # Simple Motion Definition and Simulation

using KomaMRI # hide
sys = Scanner() # hide

# This tutorial illustrates how we can add simple motion to phantoms.
# We will also see how phantoms can be stored and loaded from ``.phantom`` files.

# First, we load a static 3D brain phantom:
obj1 = brain_phantom3D()
obj1.Δw .= 0 # Removes the off-resonance
p1 = plot_phantom_map(obj1, :T2 ; height=400)
#md savefig(p1, "../assets/3-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```

# Now, we will make a copy of the phantom, to which we will add Rotation Motion
# to recreate the patient's movement inside the scanner.

#md # Note how rotations are defined with respect to the 3 axes:
#md # ```@raw html
#md # <center><img src="../../assets/head_rotation_axis.png" width="300"></center>
#md # ```

obj2 = copy(obj1)
obj2.motion = SimpleMotion([Rotation(t_start=0.0, t_end=0.25, pitch=15.0, roll=0.0, yaw=45.0)])
p2 = plot_phantom_map(obj2, :T2 ; height=400, intermediate_time_samples=4)
#md savefig(p2, "../assets/3-phantom.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```
