####################
## Phantom OBJECT ##
####################
# @everywhere begin
"""
    Phantom(x,y,ρ,T2,Δw,Dλ1,Dλ2,Dθ,ux,uy)

Phantom object.

# Arguments
 - `name`         := Name of the Phantom
 - `x`            := Spins x-coordinates.
 - `y`            := Spins y-coordinates.
 - `ρ::Matrix`    := Proton density.
 - `T2::Matrix`   := T2 map.
 - `Δw::Matrix`   := Off-resonace map;
 - `Dλ1::Matrix`  := Diffusion tensor principal eigen value map.
 - `Dλ2::Matrix`  := Diffusion tensor secondary eigen value map.
 - `Dθ::Matrix`   := Diffusion tensor angle map.
 - `ux::Function` := Displacement field x.
 - `uy::Function` := Displacement field y.
"""
struct Phantom
	name::String #Name ob the Phantom
	x::Array{Float64} #x-coordinates of spins
	y::Array{Float64} #y-coordinates of spins
	ρ::Array{Float64} #proton density
#	T1::Matrix #T1 map
	T2::Array{Float64} #T2 map
	#T2s::Array{Float64}
	Δw::Array{Float64} #Off-resonace map
	Dλ1::Array{Float64}  #Diffusion map
	Dλ2::Array{Float64}  #Diffusion map
	Dθ::Array{Float64}  #Diffusion map
	#χ::SusceptibilityModel
	ux::Function #Displacement field x
	uy::Function #Displacement field x
end
Phantom() = Phantom("spin",zeros(1,1),zeros(1,1),ones(1,1),ones(1,1),zeros(1,1),
					zeros(1,1),zeros(1,1),zeros(1,1),(x,y,t)->0,(x,y,t)->0)
size(x::Phantom) = size(x.ρ)
"""Separate object spins in a sub-group."""
Base.getindex(obj::Phantom, p::UnitRange{Int}) = begin
		Phantom(obj.name,obj.x[p],obj.y[p],obj.ρ[p],
						obj.T2[p],obj.Δw[p],obj.Dλ1[p],obj.Dλ2[p],
						obj.Dθ[p],obj.ux,obj.uy)
end
# @everywhere
Phantom(name::String,x::Array{Float64,2},y::Array{Float64,2},
	ρ::Array{Float64},T2::Array{Float64},Δw::Array{Float64},
	Dλ1::Array{Float64},Dλ2::Array{Float64},Dθ::Array{Float64},
	ux::Function,uy::Function) = begin
	x, y = x[ρ.≠0], y[ρ.≠0]
	T2 = T2[ρ.≠0]
	Δw = Δw[ρ.≠0]
	Dλ1 = Dλ1[ρ.≠0]
	Dλ2 = Dλ2[ρ.≠0]
	Dθ = Dθ[ρ.≠0]
	ρ = ρ[ρ.≠0]
	Phantom(name,x,y,ρ,T2,Δw,Dλ1,Dλ2,Dθ,ux,uy)
end
# @everywhere
Phantom(name::String,x::Array{Float64,2},y::Array{Float64,2},
	ρ::Array{Float64},T2::Array{Float64},Dλ1::Array{Float64},Dλ2::Array{Float64},
	Dθ::Array{Float64},ux::Function,uy::Function) = begin
	x, y = x[ρ.≠0], y[ρ.≠0]
	T2 = T2[ρ.≠0]
	Dλ1 = Dλ1[ρ.≠0]
	Dλ2 = Dλ2[ρ.≠0]
	Dθ = Dθ[ρ.≠0]
	ρ = ρ[ρ.≠0]
	Phantom(name,x,y,ρ,T2,zeros(size(ρ)),Dλ1,Dλ2,Dθ,ux,uy)
end
# @everywhere
Phantom(name::String,x::Array{Float64,2},y::Array{Float64,2},ρ::Array{Float64},
	T2::Array{Float64},Δw::Array{Float64}) = begin
	x, y = x[ρ.≠0], y[ρ.≠0]
	T2 = T2[ρ.≠0]
	Δw = Δw[ρ.≠0]
	ρ = ρ[ρ.≠0]
	Phantom(name,x,y,ρ,T2,Δw,zeros(size(ρ)),zeros(size(ρ)),zeros(size(ρ)),(x,y,t)->0,(x,y,t)->0)
end
+(s1::Phantom,s2::Phantom) =begin
	Phantom(s1.name*"_"*s2.name,[s1.x;s2.x],[s1.y;s2.y],[s1.ρ;s2.ρ],[s1.T2;s2.T2],
	[s1.Δw;s2.Δw],[s1.Dλ1;s2.Dλ1],[s1.Dλ2;s2.Dλ2],[s1.Dθ;s2.Dθ],s1.ux,s1.uy)
end
# Movement related commands
StartAt(s::Phantom,t0::Float64) = Phantom(s.name,s.x,s.y,s.ρ,s.T2,s.Δw,s.Dλ1,s.Dλ2,s.Dθ,(x,y,t)->s.ux(x,y,t.+t0),(x,y,t)->s.uy(x,y,t.+t0))
FreezeAt(s::Phantom,t0::Float64) = Phantom(s.name*"STILL",s.x.+s.ux(s.x,s.y,t0),s.y.+s.uy(s.x,s.y,t0),s.ρ,s.T2,s.Δw,s.Dλ1,s.Dλ2,s.Dθ,(x,y,t)->0,(x,y,t)->0)
# Getting maps
get_DxDy2D(obj::Phantom) = begin
	P(i) = rotz(obj.Dθ[i])[1:2,1:2];
	D(i) = [obj.Dλ1[i] 0;0 obj.Dλ2[i]]
	nx = [1;0]; ny = [0;1]
	Dx = [nx'*P(i)'*D(i)*P(i)*nx for i=1:prod(size(obj.Dλ1))]
	Dy = [ny'*P(i)'*D(i)*P(i)*ny for i=1:prod(size(obj.Dλ1))]
	Dx, Dy
end

"""
Heart-like LV phantom. The variable `α` is for strech, `β` for contraction, and `γ` for rotation.

 - Strech:       -α v x (FOV²/4-r²) sin(ωHR t/2)²
 - Contraction:  -β x αt
 - Rotation:      γ(x cos(θt) + y sin(θt)) 4v sin.(ωHR t)
"""
heart_phantom(x,y,FOV,α=1,β=1,γ=1,fat_bool::Bool=false) = begin
	⚪(x,y,D) =  (x.^2 .+ y.^2 .<= D^2/4)*1. #Circle of diameter D
	v = FOV/4 #m/s 1/16 th of the FOV during acquisition
	ωHR = 2π/1 #One heart-beat in one second
	θ(t) = -π/4*γ*(sin(ωHR*t).*(sin(ωHR*t).>0)+0.25.*sin(ωHR*t).*(sin(ωHR*t).<0) )
	ux(x,y,t) = begin
		-v.*(.3e5.*α.*x .*(FOV^2/4 .-(x.^2 .+y.^2)) ).*sin.(ωHR.*(t.+1/3)/2).^2 .+
		x.*v.*sin.(ωHR*t./2).^2. *β .+
		(x.*(cos.(θ.(t)).-1).+y.*sin.(θ.(t)))
	end
	uy(x,y,t) = begin
		-v.*(.3e5.*α.*y .*(FOV^2/4 .-(x.^2 .+y.^2)) ).*sin.(ωHR.*(t.+1/3)/2).^2 .+
		y.*v.*sin.(ωHR*t./2).^2. *β .+
		(y.*(cos.(θ.(t)).-1).-x.*sin.(θ.(t)))
	end
	# Water spins
	ρ = ⚪(x,y,9/10*FOV) #proton density
	# Diffusion tensor model
	D = 2e-9 #Diffusion of free water m2/s
	Dλ1 = D*⚪(x,y,FOV*6/11) .+ D*(⚪(x,y,FOV*9/10).-⚪(x,y,FOV*6/11)) #Diffusion map
	Dλ2 = D*⚪(x,y,FOV*6/11) .+ D/20*(⚪(x,y,FOV*9/10).-⚪(x,y,FOV*6/11)) #Diffusion map
	Dθ = ⚪(x,y,FOV*6/11)*0 .+ atan.(x,-y).*(⚪(x,y,FOV*9/10).-⚪(x,y,FOV*6/11)) #Diffusion map
	T1 = (1026*(⚪(x,y,9/10*FOV).-⚪(x,y,6/11*FOV)).+1400 .*⚪(x,y,6/11*FOV))/1000 #Myocardial T1
	T2 = (42*(⚪(x,y,9/10*FOV).-⚪(x,y,6/11*FOV)).+308*⚪(x,y,6/11*FOV))/1000 #T2 map [s]
	# Generating Phantoms
	heart = Phantom("SIMPLE",x,y,ρ,T2,Dλ1,Dλ2,Dθ,ux,uy)
	# Fat spins
	ρ_fat = .5*ρ.*(⚪(x,y,FOV).-⚪(x,y,9/10*FOV))
	Δw_fat = 2π*220*( ⚪(x,y,FOV).-⚪(x,y,9/10*FOV) )
	T1_fat = 800*(⚪(x,y,FOV).-⚪(x,y,9/10*FOV))/1000
	T2_fat = 120*(⚪(x,y,FOV).-⚪(x,y,9/10*FOV))/1000 #T2 map [s]
	fat = Phantom("FAT",x,y,ρ_fat,T2_fat,Δw_fat)
	#Resulting phantom
	obj = fat_bool ? heart+fat : heart #concatenating spins
end
"""
B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic digital brain phantom"
NeuroImage, in review - 2006.
B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital brain phantoms for creation of validation image data bases"
IEEE TMI, in review - 2006

https://brainweb.bic.mni.mcgill.ca/brainweb
"""
function brain_phantom2D(;axis="axial",ss=4)
    path = @__DIR__
    data = MAT.matread(path*"/exampledata/brain2D.mat")

    class = data[axis][1:ss:end,1:ss:end]
    Δx = 1e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    T2 = (class.==23)*329 .+ #CSF
        (class.==46)*83 .+ #GM
        (class.==70)*70 .+ #WM
        (class.==93)*70 .+ #FAT1
        (class.==116)*47 .+ #MUSCLE
        (class.==139)*329 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*70 .+ #FAT2
        (class.==232)*329 .+ #DURA
        (class.==255)*70 #MARROW
    T2s = (class.==23)*58 .+ #CSF
        (class.==46)*69 .+ #GM
        (class.==70)*61 .+ #WM
        (class.==93)*58 .+ #FAT1
        (class.==116)*30 .+ #MUSCLE
        (class.==139)*58 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*61 .+ #FAT2
        (class.==232)*58 .+ #DURA
        (class.==255)*61 #MARROW
        (class.==255)*70 #MARROW
    T1 = (class.==23)*2569 .+ #CSF
        (class.==46)*833 .+ #GM
        (class.==70)*500 .+ #WM
        (class.==93)*350 .+ #FAT1
        (class.==116)*900 .+ #MUSCLE
        (class.==139)*569 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*500 .+ #FAT2
        (class.==232)*2569 .+ #DURA
        (class.==255)*500 #MARROW
    ρ = (class.==23)*1 .+ #CSF
        (class.==46)*.86 .+ #GM
        (class.==70)*.77 .+ #WM
        (class.==93)*1 .+ #FAT1
        (class.==116)*1 .+ #MUSCLE
        (class.==139)*1 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*.77 .+ #FAT2
        (class.==232)*1 .+ #DURA
        (class.==255)*.77 #MARROW

    phantom = Phantom("brain2D_"*axis,x,y,ρ,T2*1e-3,zeros(size(ρ)))
end
