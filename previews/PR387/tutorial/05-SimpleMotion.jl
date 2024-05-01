using KomaMRI # hide
sys = Scanner() # hide

obj1 = brain_phantom3D()
obj1.Δw .= 0 # Removes the off-resonance
p1 = plot_phantom_map(obj1, :T2 ; height=400)
display(p1)

obj2 = copy(obj1)
obj2.motion = SimpleMotion([Rotation(t_start=0.0, t_end=0.25, pitch=15.0, roll=0.0, yaw=45.0)])
p2 = plot_phantom_map(obj2, :T2 ; height=400, intermediate_time_samples=4)
display(p2)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
