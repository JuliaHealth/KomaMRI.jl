using HDF5

function signal_jemris()
    path = @__DIR__
    sig = h5open(joinpath(path, "jemris_signals_epi_sphere_cs.h5"))["/signal/channels/00"]
    sig = sig[1,:] + 1im*sig[2,:]
    sig = sig[:]
    return sig
end

function phantom_sphere()
    path = @__DIR__
    fid = h5open(joinpath(path, "phantom_sphere.h5"), "r")
    x, y, z, ρ = fid["x"][:], fid["y"][:], fid["z"][:], fid["ρ"][:]
    T1, T2, T2s, Δw = fid["T1"][:], fid["T2"][:], fid["T2s"][:], fid["Δw"][:]
    return Phantom(; name="sphere", x, y, z, ρ, T1, T2, T2s, Δw)
end

function seq_epi_100x100_TE100_FOV230()

    # First Block (excitation)
    seq = Sequence()
    rf = RF((5.871650124959988e-5 - 0.0im)*ones(99), 9.9e-5, 5.0e-7, 0.0)
    seq += rf; seq = seq[2:end]
    seq[1].DUR[1] = 9.949999999999999e-5

    # Second Block (move to a corner of the k-space)
    gx = Grad(-0.015560481134096917, 0.0, 0.00041999999999999996, 0.00041999999999999996, 0.0)
    gy = Grad(0.013799413552738015, 0.0, 0.00037, 0.00037, 0.0)
    gz = Grad(0.0, 0.0)
    seq += Sequence(reshape([gx; gy; gz], :, 1))

    # Third Block (a delay)
    gx = Grad(0.0, 0.01517)
    seq += Sequence([gx])

    # Define blocks that will be repeated
    begin
        gx = Grad(0.010211565230481714, 0.001, 0.00028, 0.00028, 0.0)
        gy = Grad(0.0, 0.0)
        gz = Grad(0.0, 0.0)
        rf = RF(0.0, 0.0)
        adc = ADC(100, 0.00099, 0.000285, 0.0, 0.0)
        sxA = Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])
        sxB = Sequence(reshape([-gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])
        gx = Grad(0.0, 0.0)
        gy = Grad(-0.0017019276167022875, 0.0, 5.9999999999999995e-5, 5.9999999999999995e-5, 0.0)
        gz = Grad(0.0, 0.0)
        syC = Sequence(reshape([gx; gy; gz], :, 1))
    end

    # Repetition of blocks (move in the k-space)
    for _ in 1:49
        seq += sxA
        seq += syC
        seq += sxB
        seq += syC
    end

    # Final Blocks
    seq += sxA
    seq += syC
    seq += sxB

    return seq
end
