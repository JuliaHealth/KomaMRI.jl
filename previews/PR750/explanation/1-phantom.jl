using KomaMRI #hide
import DisplayAs #hide

obj = brain_phantom2D()

p1 = plot_phantom_map(obj, :T1; height=450)
display(p1);

obj.name

obj.x

obj.motion

obj[1:2000]
p2 = plot_phantom_map(obj[1:1000], :T2 ; height=450) #hide
display(p2);

obj2 = pelvis_phantom2D()
obj2.x .+= 0.1; obj.x.-= 0.1 #hide
obj_sum = obj + obj2
p3 = plot_phantom_map(obj_sum, :T1 ; height=450) #hide
display(p3);

obj_mul = 3*obj;

obj.ρ

obj_mul.ρ

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
