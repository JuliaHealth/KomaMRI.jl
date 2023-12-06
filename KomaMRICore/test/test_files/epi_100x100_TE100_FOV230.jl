function epi_100x100_TE100_FOV230()

seq = Sequence()
A, T, delay, Δf = (5.871650124959988e-5 - 0.0im)*ones(99), 9.9e-5, 5.0e-7, 0.0
rf = RF(A, T, delay, Δf)
seq += rf; seq = seq[2:end]
dur = 9.949999999999999e-5
seq[1].DUR[1] = dur

A, T, rise, fall, delay = -0.015560481134096917, 0.0, 0.00041999999999999996, 0.00041999999999999996, 0.0
gx = Grad(A, T, rise, fall, delay)
A, T, rise, fall, delay = 0.013799413552738015, 0.0, 0.00037, 0.00037, 0.0
gy = Grad(A, T, rise, fall, delay)
gz = Grad(0.0, 0.0)
seq += Sequence(reshape([gx; gy; gz], :, 1))

gx = Grad(0.0, 0.01517)
seq += Sequence([gx])

for _ in 1:49

A, T, rise, fall, delay = 0.010211565230481714, 0.001, 0.00028, 0.00028, 0.0
gx = Grad(A, T, rise, fall, delay)
gy = Grad(0.0, 0.0)
gz = Grad(0.0, 0.0)
rf = RF(0.0, 0.0)
N, T, delay, Δf, ϕ = 100, 0.00099, 0.000285, 0.0, 0.0
adc = ADC(N, T, delay, Δf, ϕ)
seq += Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])

A, T, rise, fall, delay = -0.0017019276167022875, 0.0, 5.9999999999999995e-5, 5.9999999999999995e-5, 0.0
gx = Grad(0.0, 0.0)
gy = Grad(A, T, rise, fall, delay)
gz = Grad(0.0, 0.0)
seq += Sequence(reshape([gx; gy; gz], :, 1))

A, T, rise, fall, delay = -0.010211565230481714, 0.001, 0.00028, 0.00028, 0.0
gx = Grad(A, T, rise, fall, delay)
gy = Grad(0.0, 0.0)
gz = Grad(0.0, 0.0)
rf = RF(0.0, 0.0)
N, T, delay, Δf, ϕ = 100, 0.00099, 0.000285, 0.0, 0.0
adc = ADC(N, T, delay, Δf, ϕ)
seq += Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])

A, T, rise, fall, delay = -0.0017019276167022875, 0.0, 5.9999999999999995e-5, 5.9999999999999995e-5, 0.0
gx = Grad(0.0, 0.0)
gy = Grad(A, T, rise, fall, delay)
gz = Grad(0.0, 0.0)
seq += Sequence(reshape([gx; gy; gz], :, 1))

end

A, T, rise, fall, delay = 0.010211565230481714, 0.001, 0.00028, 0.00028, 0.0
gx = Grad(A, T, rise, fall, delay)
gy = Grad(0.0, 0.0)
gz = Grad(0.0, 0.0)
rf = RF(0.0, 0.0)
N, T, delay, Δf, ϕ = 100, 0.00099, 0.000285, 0.0, 0.0
adc = ADC(N, T, delay, Δf, ϕ)
seq += Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])

A, T, rise, fall, delay = -0.0017019276167022875, 0.0, 5.9999999999999995e-5, 5.9999999999999995e-5, 0.0
gx = Grad(0.0, 0.0)
gy = Grad(A, T, rise, fall, delay)
gz = Grad(0.0, 0.0)
seq += Sequence(reshape([gx; gy; gz], :, 1))

A, T, rise, fall, delay = -0.010211565230481714, 0.001, 0.00028, 0.00028, 0.0
gx = Grad(A, T, rise, fall, delay)
gy = Grad(0.0, 0.0)
gz = Grad(0.0, 0.0)
rf = RF(0.0, 0.0)
N, T, delay, Δf, ϕ = 100, 0.00099, 0.000285, 0.0, 0.0
adc = ADC(N, T, delay, Δf, ϕ)
seq += Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])

seq += Sequence([Grad(0.0, 0.0)])
dur = 0.00011999999999999999
seq[11].DUR[1] = dur

seq += Sequence([Grad(0.0, 0.0)])
dur = 0.0015599999999999998
seq[12].DUR[1] = dur

return seq
end
