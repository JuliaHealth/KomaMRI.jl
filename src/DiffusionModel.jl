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

##Example 1
using Plots
gr()
# gererate a sample and plot it
λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x
density(s) = sqrt(λ(25,s,1)*λ(s,1,1))/s
data = getRandom(density, 4,16, 1000, 100);
histogram(data, bin=50, normalize=:pdf, xlabel="Mpipi^2 (GeV)", ylabel="Entries / 60 MeV", label="Dalitz plot projection")
x = range(4,16,length=200)
plot!(x,density.(x)/(sum(density.(x))*(x[2]-x[1])),linewidth=3)
##Example 2
using Plots
gr()
# gererate a sample and plot it
density(x) = 1-abs(x)
data = getRandom(density, -1,1, 1000, 100);
histogram(data, bin=20, normalize=:pdf, xlabel="Δx")
x = range(-1,1,length=200)
plot!(x,density.(x)/(sum(density.(x))*(x[2]-x[1])),linewidth=3)
