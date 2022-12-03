using Interpolations, Adapt, PlotlyJS

b1 = [1, 1.5, 2, 1.5, 1, 0]
dt = [1, 1, 2, 1, 1]
t = [0; cumsum(dt)]
itpp = extrapolate(interpolate((t, ), b1, Gridded(Constant{Previous}())), 0)
itp = interpolate((t, ), b1, Gridded(Constant{Previous}()))
cuitp = adapt(CuArray{Float32}, itp)
x = 0:.001:sum(dt)
y = @time itpp.(x);
y = @time cuitp.(x);
plot(x,y)
##
# using KomaMRI
# (gx,gy,gz) = get_grads(seq, t)