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
    Trf = 1e-4                     # From File: 9.9e-5
    Arf = (π/2)/(2π*γ*Trf)         # From File: 5.871650124959988e-5. 90° rect excitation pulse
    rf = RF(Arf, Trf)              # From File: RF(Arf*ones(99), Trf, 5.0e-7, 0.0)
    seq += rf;                     # From File: seq[1].DUR[1] = 9.949999999999999e-5

    # Define params
    begin
        # User params
        FOVx = 0.230                    # Set
        Nx = 100                        # Set
        FOVy = 0.230                    # Set
        Ny = 100                        # Set
        TE = 0.1                        # Set
        Δtadc = 0.00001                 # Set
        Smax = 10*1_000_000_000/2π/γ    # Set

        # For gradients which moves in the k-space (and ADC)
        Area_adc = Nx/FOVx/γ
        Tgx = Nx*Δtadc                                  # From File: 0.001
        Agx = Area_adc/Tgx                              # From File: 0.010211565230481714
        ζgx = ceil(100000*(Agx/Smax))/100000            # From File: 0.00028 (JEMRIS non-optimal definition for trapezoidal with flat-top)
        Tadc = Tgx - Δtadc                              # From File: 0.00099
        ΔDadc = ζgx + Δtadc/2                           # From File: 0.000285
        Area_gy = 1/FOVy/γ
        ζgy = ceil(100000*sqrt(Area_gy/Smax))/100000    # From File: 5.9999999999999995e-5 (JEMRIS optimal definition for trapezoidal without flat-top)
        Agy = Area_gy/ζgy                               # From File: 0.0017019276167022875

        # For first corner in the k-space (half of the covered rectangle in the k-space)
        Area_kx = Area_adc + Agx*ζgx
        Area_gxo = 0.5*Area_kx
        ζgxo = ceil(100000*sqrt(Area_gxo/Smax))/100000  # From File: 0.00041999999999999996 (JEMRIS optimal definition for trapezoidal without flat-top)
        Agxo = Area_gxo/ζgxo                            # From File: 0.015560481134096917
        Area_ky = Ny*Area_gy
        Area_gyo = 0.5*Area_ky
        ζgyo = ceil(100000*sqrt(Area_gyo/Smax))/100000  # From File: 0.00037
        Agyo = Area_gyo/ζgyo                            # From File: 0.013799413552738015 (JEMRIS optimal definition for trapezoidal without flat-top)

        # For delay (from JEMRIS center-to-center TE definition)
        ΔD = TE - (Trf/2 + 2*ζgxo + Ny*ζgx + Ny*Tgx/2 + (Ny-1)*ζgy)
    end

    # Second Block (move to a corner of the k-space)
    gx = Grad(-Agxo, 0.0, ζgxo)
    gy = Grad(Agyo, 0.0, ζgyo)
    gz = Grad(0.0, 0.0)
    seq += Sequence(reshape([gx; gy; gz], :, 1))

    # Third Block (a delay)
    gx = Grad(0.0, ΔD)
    seq += Sequence([gx])

    # Define blocks that will be repeated
    gx = Grad(Agx, Tgx, ζgx)
    gy = Grad(0.0, 0.0)
    gz = Grad(0.0, 0.0)
    rf = RF(0.0, 0.0)
    adc = ADC(Nx, Tadc, ΔDadc)
    sxA = Sequence(reshape([gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])
    sxB = Sequence(reshape([-gx; gy; gz], :, 1), reshape([rf], :, 1), [adc])
    gx = Grad(0.0, 0.0)
    gy = Grad(-Agy, 0.0, ζgy)
    gz = Grad(0.0, 0.0)
    syC = Sequence(reshape([gx; gy; gz], :, 1))

    # Repetition of blocks (move in the k-space)
    for _ in 1:(Ny÷2-1)
        seq += sxA  # moves to the right in kx and samples the signal
        seq += syC  # moves a blip down in ky
        seq += sxB  # moves to the left in kx and samples the signal
        seq += syC  # moves a blip down in ky again
    end
    seq += sxA  # moves to the right in kx and samples the signal
    seq += syC  # moves a blip down in ky
    seq += sxB  # moves to the left in kx and samples the signal

    return seq
end
