"""
    obj = Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, motion)

The Phantom struct. Most of its field names are vectors, with each element associated with
a property value representing a spin. This struct serves as an input for the simulation.

# Arguments
- `name`: (`::String`) phantom name
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `ρ`: (`::AbstractVector{T<:Real}`) spin proton density vector
- `T1`: (`::AbstractVector{T<:Real}`, `[s]`) spin T1 parameter vector
- `T2`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2 parameter vector
- `T2s`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2s parameter vector
- `Δw`: (`::AbstractVector{T<:Real}`, `[rad/s]`) spin off-resonance parameter vector
- `Dλ1`: (`::AbstractVector{T<:Real}`) spin Dλ1 (diffusion) parameter vector
- `Dλ2`: (`::AbstractVector{T<:Real}`) spin Dλ2 (diffusion) parameter vector
- `Dθ`: (`::AbstractVector{T<:Real}`) spin Dθ (diffusion) parameter vector
- `motion`: (`::AbstractMotion{T<:Real}`) motion

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = Phantom(x=[0.0])

julia> obj.ρ
```
"""
@with_kw mutable struct Phantom{T<:Real}
    name::String           = "spins"
    x                      :: AbstractVector{T}
    y::AbstractVector{T}   = zeros(eltype(x), size(x))
    z::AbstractVector{T}   = zeros(eltype(x), size(x))
    ρ::AbstractVector{T}   = ones(eltype(x), size(x))
    T1::AbstractVector{T}  = ones(eltype(x), size(x)) * 1_000_000
    T2::AbstractVector{T}  = ones(eltype(x), size(x)) * 1_000_000
    T2s::AbstractVector{T} = ones(eltype(x), size(x)) * 1_000_000
    #Off-resonance related
    Δw::AbstractVector{T} = zeros(eltype(x), size(x))
    #χ::Vector{SusceptibilityModel}
    #Diffusion
    Dλ1::AbstractVector{T} = zeros(eltype(x), size(x))
    Dλ2::AbstractVector{T} = zeros(eltype(x), size(x))
    Dθ::AbstractVector{T}  = zeros(eltype(x), size(x))
    #Diff::Vector{DiffusionModel}  #Diffusion map
    #Motion
    motion::AbstractMotion{T} = NoMotion{eltype(x)}() 
end

const NON_STRING_PHANTOM_FIELDS = Iterators.filter(x -> fieldtype(Phantom, x) != String,         fieldnames(Phantom))
const VECTOR_PHANTOM_FIELDS     = Iterators.filter(x -> fieldtype(Phantom, x) <: AbstractVector, fieldnames(Phantom))

"""Size and length of a phantom"""
size(x::Phantom) = size(x.ρ)
Base.length(x::Phantom) = length(x.ρ)
# To enable to iterate and broadcast over the Phantom
Base.iterate(x::Phantom) = (x[1], 2)
Base.iterate(x::Phantom, i::Integer) = (i <= length(x)) ? (x[i], i + 1) : nothing
Base.lastindex(x::Phantom) = length(x)
Base.getindex(x::Phantom, i::Integer) = x[i:i]
Base.view(x::Phantom, i::Integer) = @view(x[i:i])

"""Compare two phantoms"""
function Base.:(==)(obj1::Phantom, obj2::Phantom)
    if length(obj1) != length(obj2) return false end
    return reduce(&, [getfield(obj1, field) == getfield(obj2, field) for field in NON_STRING_PHANTOM_FIELDS])
end
function Base.:(≈)(obj1::Phantom, obj2::Phantom)
    if length(obj1) != length(obj2) return false end
    return reduce(&, [getfield(obj1, field)  ≈ getfield(obj2, field) for field in NON_STRING_PHANTOM_FIELDS])
end

"""Separate object spins in a sub-group"""
function Base.getindex(obj::Phantom, p)
    fields = []
    for field in NON_STRING_PHANTOM_FIELDS
        push!(fields, (field, getfield(obj, field)[p]))
    end
    return Phantom(; name=obj.name, fields...)
end
function Base.view(obj::Phantom, p)
    fields = []
    for field in NON_STRING_PHANTOM_FIELDS
        push!(fields, (field, @view(getfield(obj, field)[p])))
    end
    return Phantom(; name=obj.name, fields...)
end

"""Addition of phantoms"""
+(obj1::Phantom, obj2::Phantom) = begin
    name = first(obj1.name * "+" * obj2.name, 50) # The name is limited to 50 characters
    fields = []
    for field in VECTOR_PHANTOM_FIELDS
        push!(fields, (field, [getfield(obj1, field); getfield(obj2, field)]))
    end
    return Phantom(; 
        name = name, 
        fields..., 
        motion = vcat(obj1.motion, obj2.motion, length(obj1), length(obj2)))
end

"""Scalar multiplication of a phantom"""
*(α::Real, obj::Phantom) = begin
    obj1 = copy(obj)
    obj1.ρ .*= α
    return obj1
end

"""dims = get_dims(obj)"""
function get_dims(obj::Phantom)
    dims = Bool[]
    push!(dims, any(x -> x != 0, obj.x))
    push!(dims, any(x -> x != 0, obj.y))
    push!(dims, any(x -> x != 0, obj.z))
    return dims
end

"""
    obj = heart_phantom(
        circumferential_strain, radial_strain, rotation_angle; 
        heart_rate, asymmetry
    )

Heart-like LV 2D phantom. The variable `circumferential_strain` and `radial_strain` are for streching (if positive) 
or contraction (if negative). `rotation_angle` is for rotation.

# Keywords
- `circumferential_strain`: (`::Real`, `=-0.3`) contraction parameter. Between -1 and 1
- `radial_strain`: (`::Real`, `=-0.3`) contraction parameter. Between -1 and 1
- `rotation_angle`: (`::Real`, `=15.0`, `[º]`) maximum rotation angle
- `heart_rate`: (`::Real`, `=60`, `[bpm]`) heartbeat frequency
- `temporal_asymmetry`: (`::Real`, `=0.2`) time fraction of the period in which the systole occurs. Therefore, diastole lasts for `period * (1 - temporal_asymmetry)`

# Returns
- `obj`: (`::Phantom`) Heart-like LV phantom struct
"""
function heart_phantom(;
    circumferential_strain=-0.3,
    radial_strain=-0.3,
    rotation_angle=15.0,
    heart_rate=60,
    temporal_asymmetry=0.2,
)
    #PARAMETERS
    FOV = 10e-2 # [m] Diameter ventricule
    N = 21
    Δxr = FOV / (N - 1) #Aprox rec resolution, use Δx_pix and Δy_pix
    Ns = 50 #number of spins per voxel
    Δx = Δxr / sqrt(Ns) #spin separation
    #POSITIONS
    x = y = (-FOV / 2):Δx:(FOV / 2 - Δx) #spin coordinates
    x, y = x .+ y' * 0, x * 0 .+ y' #grid points
    #PHANTOM
    ⚪(R) = (x .^ 2 .+ y .^ 2 .<= R^2) * 1.0 #Circle of radius R
    period = 60 / heart_rate     # [s] Period
    # Water spins
    R = 9 / 10 * FOV / 2
    r = 6 / 11 * FOV / 2
    ring = ⚪(R) .- ⚪(r)
    ρ = ⚪(r) .+ 0.9 * ring #proton density
    # Diffusion tensor model
    D = 2e-9 #Diffusion of free water m2/s
    D1, D2 = D, D / 20
    Dλ1 = D1 * ⚪(R) #Diffusion map
    Dλ2 = D1 * ⚪(r) .+ D2 * ring #Diffusion map
    Dθ = atan.(x, -y) .* ring #Diffusion map
    T1 = (1400 * ⚪(r) .+ 1026 * ring) * 1e-3 #Myocardial T1
    T2 = (308 * ⚪(r) .+ 42 * ring) * 1e-3 #T2 map [s]
    # Generating Phantoms
    phantom = Phantom(;
        name="LeftVentricle",
        x=x[ρ .!= 0],
        y=y[ρ .!= 0],
        ρ=ρ[ρ .!= 0],
        T1=T1[ρ .!= 0],
        T2=T2[ρ .!= 0],
        Dλ1=Dλ1[ρ .!= 0],
        Dλ2=Dλ2[ρ .!= 0],
        Dθ=Dθ[ρ .!= 0],
        motion=MotionList(
            HeartBeat(
                circumferential_strain,
                radial_strain,
                0.0,
                Periodic(; period=period, asymmetry=temporal_asymmetry),
            ),
            Rotate(
                0.0, 0.0, rotation_angle, Periodic(; period=period, asymmetry=temporal_asymmetry)
            ),
        ),
    )
    return phantom
end

"""
    phantom = brain_phantom2D(;axis="axial", ss=4)

Creates a two-dimensional brain Phantom struct.
Default ss=4 sample spacing is 2 mm. Original file (ss=1) sample spacing is .5 mm.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb/tissue_mr_parameters.txt

# Keywords
- `axis`: (`::String`, `="axial"`, opts=[`"axial"`, `"coronal"`, `"sagittal"`]) orientation of the phantom
- `ss`: (`::Integer or ::Vector{Integer}`, `=4`) subsampling parameter for all axes if scaler, per axis if 2 element vector [ssx, ssy]
- `us`: (`::Integer or ::Vector{Integer}`, `=1`)  upsampling parameter for all axes if scaler, per axis if 2 element vector [usx, usy], if used ss is set to ss=1


# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom2D(; axis="sagittal", ss=1)

julia> obj = brain_phantom2D(; axis="axial", us=[1, 2])

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom2D(; axis="axial", ss=4, us=1)
    # check and filter input    
    ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(2, ss, us)

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path * "/phantom/brain2D.mat")

    # subsample or upsample the phantom data
    class = repeat(data[axis][1:ssx:end, 1:ssy:end]; inner=[usx, usy])

    # Define spin position vectors
    Δx = .5e-3 * ssx / usx
    Δy = .5e-3 * ssy / usy
    M, N = size(class)
    FOVx = (M - 1) * Δx #[m]
    FOVy = (N - 1) * Δy #[m]
    x = (-FOVx / 2):Δx:(FOVx / 2) #spin coordinates
    y = (-FOVy / 2):Δy:(FOVy / 2) #spin coordinates
    x, y = x .+ y' * 0, x * 0 .+ y' #grid points

    # Define spin property vectors
    T2 =
        (class .== 23) * 329 .+ #CSF
        (class .== 46) * 83 .+ #GM
        (class .== 70) * 70 .+ #WM
        (class .== 93) * 70 .+ #FAT1
        (class .== 116) * 47 .+ #MUSCLE
        (class .== 139) * 329 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 70 .+ #FAT2
        (class .== 232) * 329 .+ #DURA
        (class .== 255) * 70 #MARROW
    T2s =
        (class .== 23) * 58 .+ #CSF
        (class .== 46) * 69 .+ #GM
        (class .== 70) * 61 .+ #WM
        (class .== 93) * 58 .+ #FAT1
        (class .== 116) * 30 .+ #MUSCLE
        (class .== 139) * 58 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 61 .+ #FAT2
        (class .== 232) * 58 .+ #DURA
        (class .== 255) * 61 #MARROW
    T1 =
        (class .== 23) * 2569 .+ #CSF
        (class .== 46) * 833 .+ #GM
        (class .== 70) * 500 .+ #WM
        (class .== 93) * 350 .+ #FAT1
        (class .== 116) * 900 .+ #MUSCLE
        (class .== 139) * 569 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 500 .+ #FAT2
        (class .== 232) * 2569 .+ #DURA
        (class .== 255) * 500 #MARROW
    ρ =
        (class .== 23) * 1 .+ #CSF
        (class .== 46) * 0.86 .+ #GM
        (class .== 70) * 0.77 .+ #WM
        (class .== 93) * 1 .+ #FAT1
        (class .== 116) * 1 .+ #MUSCLE
        (class .== 139) * 1 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 0.77 .+ #FAT2
        (class .== 232) * 1 .+ #DURA
        (class .== 255) * 0.77 #MARROW
    Δw_fat = -220 * 2π
    Δw =
        (class .== 93) * Δw_fat .+ #FAT1
        (class .== 209) * Δw_fat    #FAT2
    T1 = T1 * 1e-3
    T2 = T2 * 1e-3
    T2s = T2s * 1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(;
        name="brain2D_" * axis,
        x=y[ρ .!= 0],
        y=x[ρ .!= 0],
        z=0 * x[ρ .!= 0],
        ρ=ρ[ρ .!= 0],
        T1=T1[ρ .!= 0],
        T2=T2[ρ .!= 0],
        T2s=T2s[ρ .!= 0],
        Δw=Δw[ρ .!= 0],
    )
    return obj
end

"""
    obj = brain_phantom3D(; ss=4, us=1, start_end=[160,200])

Creates a three-dimentional brain Phantom struct.
Default ss=4 sample spacing is 2 mm. Original file (ss=1) sample spacing is .5 mm. 

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb/tissue_mr_parameters.txt

# Keywords
- `ss`: (`::Integer or ::Vector{Integer}`, `=4`) subsampling parameter for all axes if scaler, per axis if 3 element vector [ssx, ssy, ssz]
- `us`: (`::Integer or ::Vector{Integer}`, `=1`)  upsampling parameter for all axes if scaler, per axis if 3 element vector [usx, usy, usz]
- `start_end`: (`::Vector{Integer}`, `=[160,200]`) z index range of presampled phantom, 180 is center

# Returns
- `obj`: (`::Phantom`) 3D Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom3D(; ss=5)

julia> obj = brain_phantom3D(; us=[2, 2, 1])

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom3D(; ss=4, us=1, start_end=[160, 200])
    # check and filter input    
    ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(3, ss, us)

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path * "/phantom/brain3D.mat")

    # subsample or upsample the phantom data
    class = repeat(
        data["data"][1:ssx:end, 1:ssy:end, start_end[1]:ssz:start_end[2]];
        inner=[usx, usy, usz],
    )

    # Define spin position vectors
    Δx = .5e-3 * ssx / usx
    Δy = .5e-3 * ssy / usy
    Δz = .5e-3 * ssz / usz
    M, N, Z = size(class)
    FOVx = (M - 1) * Δx #[m]
    FOVy = (N - 1) * Δy #[m]
    FOVz = (Z - 1) * Δz #[m]
    xx = reshape((-FOVx / 2):Δx:(FOVx / 2), M, 1, 1) #spin coordinates
    yy = reshape((-FOVy / 2):Δy:(FOVy / 2), 1, N, 1) #spin coordinates
    zz = reshape((-FOVz / 2):Δz:(FOVz / 2), 1, 1, Z) #spin coordinates
    x = 1 * xx .+ 0 * yy .+ 0 * zz
    y = 0 * xx .+ 1 * yy .+ 0 * zz
    z = 0 * xx .+ 0 * yy .+ 1 * zz

    # Define spin property vectors
    T2 =
        (class .== 23) * 329 .+ #CSF
        (class .== 46) * 83 .+ #GM
        (class .== 70) * 70 .+ #WM
        (class .== 93) * 70 .+ #FAT1
        (class .== 116) * 47 .+ #MUSCLE
        (class .== 139) * 329 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 70 .+ #FAT2
        (class .== 232) * 329 .+ #DURA
        (class .== 255) * 70 #MARROW
    T2s =
        (class .== 23) * 58 .+ #CSF
        (class .== 46) * 69 .+ #GM
        (class .== 70) * 61 .+ #WM
        (class .== 93) * 58 .+ #FAT1
        (class .== 116) * 30 .+ #MUSCLE
        (class .== 139) * 58 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 61 .+ #FAT2
        (class .== 232) * 58 .+ #DURA
        (class .== 255) * 61 #MARROW
    T1 =
        (class .== 23) * 2569 .+ #CSF
        (class .== 46) * 833 .+ #GM
        (class .== 70) * 500 .+ #WM
        (class .== 93) * 350 .+ #FAT1
        (class .== 116) * 900 .+ #MUSCLE
        (class .== 139) * 569 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 500 .+ #FAT2
        (class .== 232) * 2569 .+ #DURA
        (class .== 255) * 500 #MARROW
    ρ =
        (class .== 23) * 1 .+ #CSF
        (class .== 46) * 0.86 .+ #GM
        (class .== 70) * 0.77 .+ #WM
        (class .== 93) * 1 .+ #FAT1
        (class .== 116) * 1 .+ #MUSCLE
        (class .== 139) * 1 .+ #SKIN/MUSCLE
        (class .== 162) * 0 .+ #SKULL
        (class .== 185) * 0 .+ #VESSELS
        (class .== 209) * 0.77 .+ #FAT2
        (class .== 232) * 1 .+ #DURA
        (class .== 255) * 0.77 #MARROW
    Δw_fat = -220 * 2π
    Δw =
        (class .== 93) * Δw_fat .+ #FAT1
        (class .== 209) * Δw_fat    #FAT2
    T1 = T1 * 1e-3
    T2 = T2 * 1e-3
    T2s = T2s * 1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(;
        name="brain3D",
        x=y[ρ .!= 0],
        y=x[ρ .!= 0],
        z=z[ρ .!= 0],
        ρ=ρ[ρ .!= 0],
        T1=T1[ρ .!= 0],
        T2=T2[ρ .!= 0],
        T2s=T2s[ρ .!= 0],
        Δw=Δw[ρ .!= 0],
    )
    return obj
end

"""
    obj = pelvis_phantom2D(; ss=4, us=1)

Creates a two-dimensional pelvis Phantom struct.
Default ss=4 sample spacing is 2 mm. Original file (ss=1) sample spacing is .5 mm.

# Keywords
- `ss`: (`::Integer or ::Vector{Integer}`, `=4`) subsampling parameter for all axes if scaler, per axis if 2 element vector [ssx, ssy]
- `us`: (`::Integer or ::Vector{Integer}`, `=1`)  upsampling parameter for all axes if scaler, per axis if 2 element vector [usx, usy]

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = pelvis_phantom2D(; ss=2])

julia> obj = pelvis_phantom2D(; us=[1, 2])

julia> pelvis_phantom2D(obj, :ρ)
```
"""
function pelvis_phantom2D(; ss=4, us=1)
    # check and filter input    
    ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(2, ss, us)

    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path * "/phantom/pelvis2D.mat")

    # subsample or upsample the phantom data
    class = repeat(data["pelvis3D_slice"][1:ssx:end, 1:ssy:end]; inner=[usx, usy])

    # Define spin position vectors
    Δx = .5e-3 * ssx / usx
    Δy = .5e-3 * ssy / usy
    M, N = size(class)
    FOVx = (M - 1) * Δx             # [m]
    FOVy = (N - 1) * Δy             # [m]
    x = (-FOVx / 2):Δx:(FOVx / 2)       # spin coordinates
    y = (-FOVy / 2):Δy:(FOVy / 2)       # spin coordinates
    x, y = x .+ y' * 0, x * 0 .+ y' # grid points

    # Define spin property vectors
    ρ =
        (class .== 102) * 0.86 .+    # Fat
        (class .== 153) * 0.9 .+     # SoftTissue
        (class .== 204) * 0.4 .+     # SpongyBone
        (class .== 255) * 0.2        # CorticalBone
    T1 =
        (class .== 102) * 366 .+   # Fat
        (class .== 153) * 1200 .+   # SoftTissue
        (class .== 204) * 381 .+    # SpongyBone
        (class .== 255) * 100       # CorticalBone
    T2 =
        (class .== 102) * 70 .+    # Fat
        (class .== 153) * 80 .+     # SoftTissue
        (class .== 204) * 52 .+     # SpongyBone
        (class .== 255) * 0.3        # CorticalBone
    T2s =
        (class .== 102) * 70 .+   # Fat
        (class .== 153) * 80 .+     # SoftTissue
        (class .== 204) * 52 .+     # SpongyBone
        (class .== 255) * 0.3   # CorticalBone
    Δw_fat = -220 * 2π
    Δw = (class .== 102) * Δw_fat # FAT1
    T1 = T1 * 1e-3
    T2 = T2 * 1e-3
    T2s = T2s * 1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(;
        name="pelvis2D",
        x=y[ρ .!= 0],
        y=x[ρ .!= 0],
        z=0 * x[ρ .!= 0],
        ρ=ρ[ρ .!= 0],
        T1=T1[ρ .!= 0],
        T2=T2[ρ .!= 0],
        T2s=T2s[ρ .!= 0],
        Δw=Δw[ρ .!= 0],
    )
    return obj
end

"""
    ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(nd, ss, us)

Utility function to check the arguments of phantom generating functions.
# Arguments
- `nd` : (`::Integer`) dimensionality of the phantom
- `ss` : (`::Integer or ::Vector{Integer}`) subsampling parameter for all axes if scaler, per axis if a 2 or 3 element vector
- `us` : (`::Integer or ::Vector{Integer}`)  upsampling parameter for all axes if scaler, per axis if a 2 or 3 element vector

# Returns
- `ssx, ssy, ssz`: (`::Integer`) valid subsampling parameters per axis
- `usx, usy, usz`: (`::Integer`)  valid upsampling parameters per axis

# Examples
```julia-repl
julia> ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(2, 1, 1)

julia> ssx, ssy, ssz, usx, usy, usz = check_phantom_arguments(3, 4, [2, 2, 2])
```
"""
function check_phantom_arguments(nd, ss, us)
    # check for valid input    
    ssz = -9999
    usz = -9999
    if length(us) > 1 || prod(us) > 1
        @info "setting ss=1 since us=$(us) defined"
        ss = 1
    end
    if nd == 3
        @assert length(ss) <= 3 "ss=$(ss) invalid, ss can have up to three components [ssx, ssy, ssz] for a 3D phantom"
        @assert length(us) <= 3 "us=$(us) invalid, us can have up to three components [usx, usy, usz] for a 3D phantom"
        if length(us) == 1
            usx = us[1]
            usy = us[1]
            usz = us[1]
        elseif length(us) == 2
            usx = us[1]
            usy = us[2]
            usz = us[2]
            @warn "Using us=$([usx, usy, usz]) in place of us=$([usx, usy])."
        else
            usx = us[1]
            usy = us[2]
            usz = us[3]
        end
        if length(ss) == 1
            ssx = ss[1]
            ssy = ss[1]
            ssz = ss[1]
        elseif length(ss) == 2
            ssx = ss[1]
            ssy = ss[2]
            ssz = ss[2]
            @warn "Using ss=$([ssx, ssy, ssz]) in place of ss=$([ssx, ssy])."
        else
            ssx = ss[1]
            ssy = ss[2]
            ssz = ss[3]
        end
    elseif nd == 2
        @assert length(ss) <= 2 "ss=$(ss) invalid, ss can have up to two components [ssx, ssy] for a 2D phantom"
        @assert length(us) <= 2 "us=$(us) invalid, us can have up to two components [usx, usy] for a 2D phantom"
        if length(us) == 1
            usx = us[1]
            usy = us[1]
        else
            usx = us[1]
            usy = us[2]
        end
        if length(ss) == 1
            ssx = ss[1]
            ssy = ss[1]
        else
            ssx = ss[1]
            ssy = ss[2]
        end
    end
    return ssx, ssy, ssz, usx, usy, usz
end
