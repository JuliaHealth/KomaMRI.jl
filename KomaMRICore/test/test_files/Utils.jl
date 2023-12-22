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

    # Gyromagnetic constant
    γ_rad_us_mT = 2π*γ*1e-9        # γ in [Hz/T]. γ_rad_us_mT in [rad/us/mT]

    # For excitation
    Trf = 100e-6                   # RF duration 100us
    Arf = (π/2)/(2π*γ*Trf)         # RF Amplitude for 90° rect excitation pulse

    # User params (defined in the original JEMRIS XML file)
    FOVx = 230e-3                       # 230mm Field-of-View in x
    Nx = 100                            # Number of points in x
    FOVy = 230e-3                       # 230mm Field-of-View in y
    Ny = 100                            # Number of points in y
    TE = 100e-3                         # 100ms Time to Echo
    Δtadc = 10e-6                       # 10us ADC sampling time
    Smax_jemris = 10                    # Gradient slew rate in [rad/m/us/ms] (These are the peculiar units used in JEMRIS. Note that the gradient units in JEMRIS are measured in [rad/m/us])
    Smax = Smax_jemris / γ_rad_us_mT    # Gradient slew rate in [T/m/s]
    Δtgr_digits = 5                     # Time resolution for gradients defined in JEMRIS source code (10us)

    # For EPI
    FOVkx = Nx / FOVx                                       # The FOV of the k-space traversed by the ADC in the x-direction
    Area_adc = FOVkx / γ                                    # The area of the rectangular gradient traversed by the ADC in the x-direction
    Tgx = Nx*Δtadc                                          # The duration of the flat top of the x-gradient
    Agx = Area_adc / Tgx                                    # The amplitude of the x-gradient
    ζgx = ceil(Agx / Smax; digits=Δtgr_digits)              # The rise and fall time of the the x-gradient
    Tadc = Tgx - Δtadc                                      # The duration of the ADC
    ΔDadc = ζgx + Δtadc/2                                   # The delay of the ADC
    Δky = 1 / FOVy                                          # The resolution of the k-space in the y-direction
    Area_gy = Δky / γ                                       # The area of the triangular blip gradient traversed in the y-direction
    ζgy = ceil(sqrt(Area_gy / Smax); digits=Δtgr_digits)    # The rise and fall time of the the triangular y-gradient (ζ = sqrt(Area/Slewrate) = sqrt((G*ζ)/(G/ζ)))
    Agy = Area_gy / ζgy                                     # The amplitude of the triangular y-gradient

    # For dephasers
    Area_gx = Area_adc + Agx * ζgx                              # The total area of the x-gradient (Area_gx = ◢ + Area_adc + ◣ = Area_adc + ■ = Area_adc + Agx * ζgx)
    Area_gxo = 1/2 * Area_gx                                    # The half-area of the x-gradient
    ζgxo = ceil(sqrt(Area_gxo / Smax); digits=Δtgr_digits)      # The rise and fall time of the the triangular dephasing x-gradient (ζ = sqrt(Area/Slewrate) = sqrt((G*ζ)/(G/ζ)))
    Agxo = Area_gxo / ζgxo                                      # The amplitude of the triangular dephasing x-gradient
    Area_ky = Ny*Area_gy                                        # The width of the k-space traversed by the y-gradient in the y-direction
    Area_gyo = 1/2 * Area_ky                                    # The half-width of the k-space in the y-direction
    ζgyo = ceil(sqrt(Area_gyo / Smax); digits=Δtgr_digits)      # The rise and fall time of the the triangular dephasing y-gradient (ζ = sqrt(Area/Slewrate) = sqrt((G*ζ)/(G/ζ)))
    Agyo = Area_gyo / ζgyo                                      # The amplitude of the triangular dephasing y-gradient

    # Define the sequence
    seq_excitation = Sequence([Grad(0.0, 0.0)], [RF(Arf, Trf)])
    seq_dephaser = Sequence(reshape([Grad(-Agxo, 0.0, ζgxo); Grad(Agyo, 0.0, ζgyo)], :, 1))
    seq_epi = Sequence(vcat(
        [(i % 2 == 0) ? Grad(Agx*(-1)^(i÷2), Tgx, ζgx) : Grad(0.0, 0.0, ζgy) for i in 0:2*Ny-1],
        [(i % 2 == 1) ? Grad(-Agy, 0.0, ζgy)           : Grad(0.0, Tgx, ζgx) for i in 0:2*Ny-1]))
    seq_epi.ADC = [(i % 2 == 0) ? ADC(Nx, Tadc, ΔDadc) : ADC(0, 0.0) for i in 0:2*Ny-1]
    seq_epi = seq_epi[1:end-1]          # Remove unnecessary last blip
    delay = Delay(TE - 1/2 * dur(seq_excitation) - dur(seq_dephaser) - 1/2 * dur(seq_epi))

    # Return the sequence
    seq = seq_excitation + seq_dephaser + delay + seq_epi
    return seq
end
