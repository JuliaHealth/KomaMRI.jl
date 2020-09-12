# UNDER DEVELOPMENT!
using StatsBase

function invertLinear(t,xmin,xmax,ymin,ymax)
    # special case, if Pol1
    ymin == ymax && return xmin + t*(xmax-xmin)
    # otherwise, Pol2
    A = (ymax-ymin)/(xmax-xmin)/2;
    B = ymin;
    C = - t*(ymin*(xmax-xmin)+(xmax-xmin)*(ymax-ymin)/2);
    x = xmin + (-B + sqrt(B^2-4*A*C))/(2*A)
    return x
end

function getRandom(func, xmin, xmax, Nvalues, Ndiv)
    x = range(xmin,xmax,length=Ndiv)
    P = func.(x)
    # remove interval where the function is negative
    P = [v<0 ? 0 : v for v in P]
    # weight is given by integral
    Pa = [(P[i]+P[i+1])/2 for i in 1:Ndiv-1]
    Pa = Pa/sum(Pa)
    # generate set on random bin indexes
    inds = StatsBase.sample(1:Ndiv-1, Weights(Pa), Nvalues)
    # convert the set of indexes to random variables inside a bin
    return [invertLinear(rand(),x[i],x[i+1],P[i],P[i+1]) for i in inds]
end

# ##Example 1
# using Plots
# gr()
# # gererate a sample and plot it
# λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x
# density(s) = sqrt(λ(25,s,1)*λ(s,1,1))/s
# data = getRandom(density, 4,16, 1000, 100);
# histogram(data, bin=50, normalize=:pdf, xlabel="Mpipi^2 (GeV)", ylabel="Entries / 60 MeV", label="Dalitz plot projection")
# x = range(4,16,length=200)
# plot!(x,density.(x)/(sum(density.(x))*(x[2]-x[1])),linewidth=3)
# ##Example 2
# using Plots
# gr()
# # gererate a sample and plot it
# density(x) = 1-abs(x)
# data = getRandom(density, -1,1, 1000, 100);
# histogram(data, bin=20, normalize=:pdf, xlabel="Δx")
# x = range(-1,1,length=200)
# plot!(x,density.(x)/(sum(density.(x))*(x[2]-x[1])),linewidth=3)

##
using LinearAlgebra, Plots
function Planes(M=10)
    #Normalization
    ϵ(m) = m==0 ? 1 : sqrt(2)
    #Eigen values
    λ(m) = π^2*m^2
    Λ =  diagm([λ(m) for m=0:M])
    #Restriction
    B = [m!=n ? ((-1)^(m+n)-1)*ϵ(m)*ϵ(n)*(λ(m)+λ(n))/(λ(m)-λ(n))^2 : 1/2 for m=0:M, n=0:M]
    #Permeability
    Bs = [ϵ(m)*ϵ(n)*(1+(-1)^(m+n)) for m=0:M, n=0:M]
    #Initial conditions and coil sensitivities in eigen basis
    U = [m==0 ? 1 : 0 for m=0:M]
    Λ, B, Bs, U 
end
function Disk(M=10)
    #J'n(αnk) = 0
    α = [0 1.841184 3.054237 3.831706 4.201189 5.317553 5.331443 6.415616 6.706133 7.015587
        7.501266 8.015237 8.536316 8.577836 9.282396 9.647422 9.969468 10.17347 10.51986 
        10.71143 11.34592]
    #Eigen values
    λ(m) = α[m+1]^2
    Λ =  diagm([λ(m) for m=0:M])
    #Restriction
    B = [m!=n ? ((-1)^(m+n)-1)*ϵ(m)*ϵ(n)*(λ(m)+λ(n))/(λ(m)-λ(n))^2 : 1/2 for m=0:M, n=0:M]
    #Permeability
    Bs = [ϵ(m)*ϵ(n)*(1+(-1)^(m+n)) for m=0:M, n=0:M]
    #Initial conditions and coil sensitivities in eigen basis
    U = [m==0 ? 1 : 0 for m=0:M]
    Λ, B, Bs, U 
end

Λ, B, _, U = Planes(100)
𝒊 = 1im;
T = 60; δ = T/2;

# r: dimensionless restriction coefficient
Eq(p,q) = U'*exp(-(p*Λ .+ 𝒊*q*B)*δ/T)*exp(-p*Λ*(T-2δ)/T)*exp(-(p*Λ .- 𝒊*q*B)*δ/T)*U
Eqasy(p,q) = (p/q)^(1/3)*exp(-1.0188/2*(p*q^2)^(1/3))
E = [Eq(.1,q) for q = 1:100]
Ea = [Eqasy(.1,q) for q = 1:100]
plot(abs.(E),label="r=0.1",yaxis=:log)
plot!(abs.(Ea),label="r=0.1",yaxis=:log)
E = [Eq(1,q) for q = 1:100]
Ea = [Eqasy(1,q) for q = 1:100]
plot!(abs.(E),label="r=1",yaxis=:log)
plot!(abs.(Ea),label="r=1",yaxis=:log)
E = [Eq(10,q) for q = 1:100]
Ea = [Eqasy(10,q) for q = 1:100]
plot!(abs.(E),label="r=10",yaxis=:log)
plot!(abs.(Ea),label="r=10",yaxis=:log)
ylims!(1e-5,1)
