using MRIsim, JLD2, MAT
fantoma = matread("C:/Users/ccp/Downloads/fantoma.mat")["fantoma"]

ρ = fantoma["PD"]
Δx = 2e-3
M, N = size(ρ)
FOVx = (M-1)*Δx #[m]
FOVy = (N-1)*Δx #[m]
xx = -FOVx/2:Δx:FOVx/2 #spin coordinates
yy = -FOVy/2:Δx:FOVy/2 #spin coordinates
x = xx   .+ yy'*0 
y = xx*0 .+ yy' #grid points
phantom = Phantom(name="brain",
                  x=y[ρ.!=0],
                  y=-x[ρ.!=0],
                  z=z[ρ.!=0],
                  ρ=fantoma["PD"][ρ.!=0], 
                  T1=fantoma["T1"][ρ.!=0]*1e-3,
                  T2=fantoma["T2"][ρ.!=0]*1e-3,
                  Δw=2π*fantoma["df"][ρ.!=0]
                  )
@save "C:/Users/ccp/Downloads/brain.phantom" phantom #Hacer ]up para version MRIsim v0.3.6 