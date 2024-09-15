using HDF5

function signal_sphere_jemris()
    path = @__DIR__
    sig = h5open(joinpath(path, "jemris_signals_epi_sphere_cs.h5"))["/signal/channels/00"]
    sig = sig[1,:] + 1im*sig[2,:]
    sig = sig[:]
    return sig
end

function signal_brain_motion_jemris()
    path = @__DIR__
    sig = h5open(joinpath(path, "jemris_signals_epi_brain_motion.h5"))["/signal/channels/00"]
    sig = sig[1,:] + 1im*sig[2,:]
    sig = sig[:]
    return sig
end

function phantom_brain_simple_motion()
    obj = phantom_brain()
    obj.motion = MotionList(Translate(0.0, 1.0, 0.0, TimeRange(0.0, 10.0)))
    return obj
end

function phantom_brain_arbitrary_motion()
    obj = phantom_brain()
    Ns = length(obj)
    t_start = 0.0
    t_end = 10.0
    dx = zeros(Ns, 2)  
    dz = zeros(Ns, 2)  
    dy = [zeros(Ns,1) ones(Ns,1)]
    obj.motion = MotionList(Path(
        dx,
        dy,
        dz,
        TimeRange(t_start, t_end)))
    return obj
end

function phantom_sphere()
    path = @__DIR__
    fid = h5open(joinpath(path, "phantom_sphere.h5"), "r")
    x, y, z, ρ = fid["x"][:], fid["y"][:], fid["z"][:], fid["ρ"][:]
    T1, T2, T2s, Δw = fid["T1"][:], fid["T2"][:], fid["T2s"][:], fid["Δw"][:]
    return Phantom(; name="sphere", x, y, z, ρ, T1, T2, T2s, Δw)
end

function phantom_brain()
    path = @__DIR__
    fid = h5open(joinpath(path, "../../../examples/2.phantoms/brain.h5"), "r")
    data = read(fid["sample/data"])
    Δx = read(fid["sample/resolution"]) * 1e-3 #[m]
    offset = read(fid["sample/offset"]) * 1e-3 #[m]
    mask = data[1, :, :, :] .!= 0
    #Maps
    ρ = data[1, :, :, :]
    T1 = 1e-3 ./ data[2, :, :, :]
    T2 = 1e-3 ./ data[3, :, :, :]
    T2s = 1e-3 ./ data[4, :, :, :]
    Δw = data[5, :, :, :]
    #Positions
    X, Y, Z = size(ρ)
    FOVx = (X - 1) * Δx[1] #[m]
    FOVy = (Y - 1) * Δx[2] #[m]
    FOVz = (Z - 1) * Δx[3] #[m]
    xx = reshape(((-FOVx / 2):Δx[1]:(FOVx / 2)), :, 1, 1) #[(-FOVx/2:Δx[1]:FOVx/2)...;]
    yy = reshape(((-FOVy / 2):Δx[2]:(FOVy / 2)), 1, :, 1) #[(-FOVy/2:Δx[2]:FOVy/2)...;;]
    zz = reshape(((-FOVz / 2):Δx[3]:(FOVz / 2)), 1, 1, :) #[(-FOVz/2:Δx[3]:FOVz/2)...;;;]
    x = xx * 1 .+ yy * 0 .+ zz * 0 .+ offset[1] #spin x coordinates
    y = xx * 0 .+ yy * 1 .+ zz * 0 .+ offset[2] #spin y coordinates
    z = xx * 0 .+ yy * 0 .+ zz * 1 .+ offset[3] #spin z coordinates
    v = 0 # m/s

    obj = Phantom(;
        name="brain",
        x=x[mask],
        y=y[mask],
        z=z[mask],
        ρ=ρ[mask],
        T1=T1[mask],
        T2=T2[mask],
        T2s=T2s[mask],
        Δw=Δw[mask],
    )

    return obj
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

    # Excitation
    ex = Sequence([Grad(0.0, 0.0)], [RF(Arf, Trf)])
    # Dephasing gradients (moving to the corner of k-space)
    dephaser = Sequence([Grad(-Agxo, 0.0, ζgxo); Grad(Agyo, 0.0, ζgyo);;])
    # Readout
    gx = Grad(Agx, Tgx, ζgx)
    adc = ADC(Nx, Tadc, ΔDadc)
    readout = Sequence([gx;;], [RF(0.0, 0.0);;], [adc])
    # Blip
    gy = -Grad(Agy, 0.0, ζgy)
    blip_neg = Sequence([Grad(0.0, 0.0); gy;;])
    # Sequence generation
    epi = Sequence()
    for i in 0:Ny-1
        epi += readout * (-1)^i
        epi += blip_neg
    end
    epi = epi[1:end-1]          # Remove unnecessary last blip
    # Calculating delayTE
    delayTE = Delay(TE - 1/2 * dur(ex) - dur(dephaser) - 1/2 * dur(epi))

    # Return the sequence
    seq = ex + dephaser + delayTE + epi
    return seq
end

function NRMSE(x, x_true) 
    return sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
end