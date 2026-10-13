const BIFTI_MM = 1e-3            # BIfTI affines are in mm
const BIFTI_ADC_UNIT = 1e-9      # BIfTI ADC unit, 10⁻³ mm²/s, in m²/s
const NO_RELAXATION = 1e6        # [s] Phantom default for T1, T2 and T2s
const GYRO_RTOL = 1e-3           # γ mismatch tolerated as the same nucleus (¹H in water vs. free)

"""
    obj = read_phantom_bifti(filename; density_threshold=0.0)
    obj = read_phantom_bifti(phantom::BiftiPhantoms.VoxelPhantom; name="BIfTI", density_threshold=0.0)

Returns the Phantom struct from a [BIfTI](https://github.com/mrx-org/bifti-phantoms)
phantom: a `.json` file defining tissues, referencing NIfTI files for the per-voxel data.
`reslice_to` and `func` mappings are applied by BiftiPhantoms.jl while loading.

Every tissue becomes one spin per voxel with `density > density_threshold`, at the voxel
centre in scanner coordinates (the phantom's `patient` position applied), with
`ρ = density`, `T1`, `T2`, `T2s = 1 / (1/T2 + 1/T2')`, `Δw = 2π·dB0` and
`Dλ1 = Dλ2 = ADC`. Tissues overlapping in a voxel give overlapping spins. Missing
relaxation (`Inf`) is stored as the `Phantom` default.

Koma has no transmit-field model, so `B1+` is ignored with a warning. `B1-` is turned
into a receiver by [`read_coil_sens_bifti`](@ref).

# Arguments
- `filename`: (`::String`) the absolute or relative path of the phantom file `.json`

# Keywords
- `density_threshold`: (`::Real`, `=0.0`) voxels with a tissue density at or below this are dropped

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = read_phantom_bifti("subj42-3T.json")

julia> plot_phantom_map(obj, :T1)
```
"""
function read_phantom_bifti(filename::AbstractString; density_threshold=0.0)
    return read_phantom_bifti(load_bifti(filename); name=basename(filename), density_threshold)
end

function read_phantom_bifti(phantom::VoxelPhantom; name="BIfTI", density_threshold=0.0)
    warn_unsupported_bifti(phantom)
    spins = [
        bifti_tissue_spins(tissue, scanner_affine(phantom, tissue_name), density_threshold)
        for (tissue_name, tissue) in phantom.tissues
    ]
    field(key) = reduce(vcat, getindex.(spins, key))
    return Phantom(;
        name,
        x=field(:x), y=field(:y), z=field(:z),
        ρ=field(:ρ), T1=field(:T1), T2=field(:T2), T2s=field(:T2s), Δw=field(:Δw),
        Dλ1=field(:D), Dλ2=field(:D), Dθ=zero(field(:D)),
    )
end

function bifti_tissue_spins(tissue, affine, density_threshold)
    voxels = findall(>(density_threshold), tissue.density)
    x, y, z = bifti_voxel_positions(affine, voxels)
    relaxation(T) = min.(T[voxels], NO_RELAXATION)
    return (;
        x, y, z,
        ρ=tissue.density[voxels],
        T1=relaxation(tissue.T1),
        T2=relaxation(tissue.T2),
        T2s=relaxation(@. 1 / (1 / tissue.T2 + 1 / tissue.T2dash)),
        Δw=2π .* tissue.dB0[voxels],
        D=BIFTI_ADC_UNIT .* tissue.ADC[voxels],
    )
end

# Voxel centres in scanner coordinates [m]. Spins and the coil grid of
# `read_coil_sens_bifti` share this expression, so spins land exactly on grid nodes.
bifti_voxel_positions(affine, voxels) = ntuple(
    axis -> [BIFTI_MM * (sum(affine[axis, d] * (I[d] - 1) for d in 1:3) + affine[axis, 4]) for I in voxels],
    3,
)

function warn_unsupported_bifti(phantom::VoxelPhantom)
    is_uniform(channels) = length(channels) == 1 && all(isone, only(channels))
    any(!is_uniform(tissue.B1_tx) for tissue in values(phantom.tissues)) &&
        @warn "BIfTI: B1+ is ignored, Koma has no transmit-field model."
    any(!is_uniform(tissue.B1_rx) for tissue in values(phantom.tissues)) &&
        @info "BIfTI: B1- is not part of the Phantom, use `read_coil_sens_bifti` for the receiver."
    gyro = phantom.config.system.gyro * 1e6
    isapprox(gyro, γ; rtol=GYRO_RTOL) ||
        @warn "BIfTI: the phantom's gyro = $gyro Hz/T is ignored, Koma simulates with γ = $γ Hz/T."
    return nothing
end

"""
    receiver = read_coil_sens_bifti(filename)
    receiver = read_coil_sens_bifti(phantom::BiftiPhantoms.VoxelPhantom)

Returns the `B1-` receive sensitivities of a [BIfTI](https://github.com/mrx-org/bifti-phantoms)
phantom as an `ArbitraryCoilSens`, one coil per `B1-` channel, on the phantom's voxel grid
in scanner coordinates. Where tissues overlap, their maps are averaged weighted by
density; outside any tissue, unweighted.

All tissues must share one axis-aligned grid (use `reslice_to` otherwise) and the same
number of `B1-` channels.

# Arguments
- `filename`: (`::String`) the absolute or relative path of the phantom file `.json`

# Returns
- `receiver`: (`::ArbitraryCoilSens`) receive coil sensitivities

# Examples
```julia-repl
julia> obj = read_phantom_bifti("phantom.json");

julia> sys = Scanner(; receiver=read_coil_sens_bifti("phantom.json"));
```
"""
read_coil_sens_bifti(filename::AbstractString) = read_coil_sens_bifti(load_bifti(filename))

function read_coil_sens_bifti(phantom::VoxelPhantom)
    tissues = collect(values(phantom.tissues))
    grid = first(tissues)
    all(t -> size(t) == size(grid) && t.affine == grid.affine, tissues) ||
        throw(ArgumentError("BIfTI: all tissues must share one grid for B1-, set `reslice_to`"))
    n_coils = length(grid.B1_rx)
    all(t -> length(t.B1_rx) == n_coils, tissues) ||
        throw(ArgumentError("BIfTI: all tissues must have the same number of B1- channels"))

    total = sum(t.density for t in tissues)
    uniform_weight = 1 / length(tissues)
    weights = [@. ifelse(total > 0, t.density / total, uniform_weight) for t in tissues]
    sens = stack([sum(w .* t.B1_rx[coil] for (w, t) in zip(weights, tissues)) for coil in 1:n_coils])

    # Reorder the voxel grid into increasing scanner x, y, z, as ArbitraryCoilSens expects.
    affine = scanner_affine(phantom, first(keys(phantom.tissues)))
    voxel_axis = [findall(!iszero, affine[axis, 1:3]) for axis in 1:3]
    all(==(1) ∘ length, voxel_axis) && allunique(only.(voxel_axis)) ||
        throw(ArgumentError("BIfTI: B1- needs a grid aligned with the scanner axes, set `reslice_to`"))
    sens = permutedims(sens, (only.(voxel_axis)..., 4))
    coordinates = Vector{Vector{Float64}}(undef, 3)
    for axis in 1:3
        d = only(voxel_axis[axis])
        line = [CartesianIndex(ntuple(k -> k == d ? i : 1, 3)) for i in 1:size(grid)[d]]
        coordinates[axis] = bifti_voxel_positions(affine, line)[axis]
        if affine[axis, d] < 0
            reverse!(coordinates[axis])
            sens = reverse(sens; dims=axis)
        end
    end
    return ArbitraryCoilSens(coordinates..., complex(sens))
end
