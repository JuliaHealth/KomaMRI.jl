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
    seq += rf; seq = seq[2:end]    # From File: seq[1].DUR[1] = 9.949999999999999e-5

    # Define params
    begin
        # Image params
        FOVx = 0.230                    # Set
        Nx = 100                        # Set
        FOVy = 0.230                    # Set
        Ny = 100                        # Set

        # For gradients which moves in the k-space (and ADC)
        Tgx = 0.001                     # Set
        ζgx = 0.00028                   # Set
        Δtadc_min = 0.00001             # Set
        Agx = Nx/(γ*Tgx*FOVx)           # From File: 0.010211565230481714
        Tadc = Tgx - Δtadc_min          # From File: 0.00099
        ΔDadc = ζgx + Δtadc_min/2       # From File: 0.000285
        ζgy = 6e-5                      # From File: 5.9999999999999995e-5
        Agy = 1/(γ*ζgy*FOVy)            # From File: 0.0017019276167022875

        # For first corner in the k-space
        ζgyo = 0.00037                  # Set
        ζgxo = 0.00042                  # From File: 0.00041999999999999996
        Agyo = 0.5*Ny/(γ*ζgyo*FOVy)     # From File: 0.013799413552738015
        Agxo = 0.5*Agx*(Tgx + ζgx)/ζgxo # From File: 0.015560481134096917
    end

    # Second Block (move to a corner of the k-space)
    gx = Grad(-Agxo, 0.0, ζgxo)
    gy = Grad(Agyo, 0.0, ζgyo)
    gz = Grad(0.0, 0.0)
    seq += Sequence(reshape([gx; gy; gz], :, 1))

    # Third Block (a delay)
    gx = Grad(0.0, 0.01517)
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
    for _ in 1:(round(Int, Ny/2)-1)
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
