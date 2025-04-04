using KomaMRI # hide
obj = brain_phantom2D() # hide

obj.motion = NoMotion()

obj.motion = Motion(Translate(0.0, 0.1, 0.2), TimeRange(0.0, 1.0), AllSpins())

obj.motion = MotionList(
    Motion(Translate(0.0, 0.1, 0.2), TimeRange(0.0, 1.0), AllSpins()),
    Motion(Rotate(0.0, 0.0, 45.0), Periodic(1.0, 0.5), SpinRange(1:1000))
)

L = 1e-3 # hide
Δx = 20e-6 # hide
x = y = z = -L/2:Δx:L/2 # hide
xx = reshape(x, (length(x),1,1)) # hide
yy = reshape(y, (1,length(y),1)) # hide
zz = reshape(z, (1,1,length(z))) # hide
x = 1*xx .+ 0*yy .+ 0*zz # hide
y = 0*xx .+ 1*yy .+ 0*zz # hide
z = 0*xx .+ 0*yy .+ 1*zz # hide
◼(L) =((abs.(x) .<= L/2) .& (abs.(y) .<= L/2) .& (abs.(z) .<= L/2)) # hide
cube = ◼(L) - ◼(L - Δx) # Hollow cube # hide
ρ = 1.0*cube #proton density # hide
T1 = copy(ρ) # hide
T1s = [100, 500, 1000, 2500, 2000, 1500] .* 1e-3 # hide
idx_T1 = 1 # hide
ϵ = 1e-5 # hide
for (i, x) in enumerate([x,y,z]) # hide
    for (j, L) in enumerate([-L/2, L/2]) # hide
        T1[(L - ϵ) .<= x .<= (L + ϵ)] .= T1s[idx_T1] # hide
        global idx_T1 += 1 # hide
    end # hide
end # hide
obj = Phantom( x=x[ρ .!= 0], y=y[ρ .!= 0], z=z[ρ .!= 0], T1 = T1[ρ .!= 0] ) # hide

p = plot_phantom_map(obj, :T1; height=450) #hide

display(p)

obj.motion = Translate(-5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), AllSpins())

p1 = plot_phantom_map(obj, :T1; time_samples=11)

display(p1)

obj.motion = Rotate(0.0, 90.0, 75.0, TimeRange(0.0, 1.0), AllSpins())

p2 = plot_phantom_map(obj, :T1; time_samples=11) #hide

display(p2)

obj.motion = Translate(-5e-4, 6e-4, 7e-4, TimeRange(0.0, 1.0), SpinRange(7500:15002))

p3 = plot_phantom_map(obj, :T1; time_samples=11) #hide

display(p3)

obj1 = brain_phantom2D() #hide
obj2 = copy(obj1) #hide
obj1.x .-= 20e-2; obj2.x .-= 20e-2 #hide
obj1.y .+= 12e-2; obj2.y .-= 12e-2 #hide

obj1.motion = MotionList(
    Translate(40e-2, 0.0, 0.0, TimeRange(0.0, 0.5),AllSpins()),
    Rotate(0.0, 0.0, 90.0, TimeRange(0.5, 1.0),AllSpins()),
)

obj2.motion = MotionList(
    Translate(40e-2, 0.0, 0.0, TimeRange(0.0, 1.0),AllSpins()),
    Rotate(0.0, 0.0, 90.0, TimeRange(0.0, 1.0),AllSpins()),
)

obj = obj1 + obj2

p4 = plot_phantom_map(obj, :T1; time_samples=11, view_2d=true) # hide

display(p4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
