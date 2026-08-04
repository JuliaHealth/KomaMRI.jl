# GPU-only receive-sensitivity support. Device-side sections run inside kernels;
# the remaining sections prepare or materialize their inputs before kernel launch.
@adapt_structure ArbitraryCoilSens

# GPU receiver representations are constructed before launching a simulation kernel.
struct BirdcageCoilGeometry{X<:AbstractVector,Y<:AbstractVector,T<:Real}
    x::X
    y::Y
    L::T
end
@adapt_structure BirdcageCoilGeometry

function gpu_birdcage_geometry(receiver, T, backend)
    ncoils = get_n_coils(receiver)
    ϕ = T(2π) .* T.(0:(ncoils - 1)) ./ T(ncoils)
    radius = T(receiver.radius)
    return BirdcageCoilGeometry(
        gpu(radius .* cos.(ϕ), backend),
        gpu(radius .* sin.(ϕ), backend),
        T(receiver.L),
    )
end

# Device-side interpolation used only while evaluating arbitrary maps in GPU kernels.
@inline function gpu_grid_interval(grid, query)
    n = length(grid)
    @inbounds begin
        if query < grid[1] || query > grid[n]
            return 1, zero(query), false
        elseif n == 1
            return 1, zero(query), query == grid[1]
        end

        low = 1
        high = n
        while high - low > 1
            middle = (low + high) >> 1
            if grid[middle] <= query
                low = middle
            else
                high = middle
            end
        end
        weight = (query - grid[low]) / (grid[high] - grid[low])
        return low, weight, true
    end
end

@inline gpu_lerp(a, b, weight) = a + (b - a) * weight

@inline function gpu_interpolate_sensitivity(
    x, y, z, coil, x_grid, y_grid, z_grid, coil_sens,
)
    ix, wx, inside_x = gpu_grid_interval(x_grid, x)
    iy, wy, inside_y = gpu_grid_interval(y_grid, y)
    iz, wz, inside_z = gpu_grid_interval(z_grid, z)
    inside_x && inside_y && inside_z || return zero(eltype(coil_sens))
    ix1 = min(ix + 1, length(x_grid))
    iy1 = min(iy + 1, length(y_grid))
    iz1 = min(iz + 1, length(z_grid))
    @inbounds begin
        c00 = gpu_lerp(
            coil_sens[ix, iy, iz, coil], coil_sens[ix1, iy, iz, coil], wx,
        )
        c10 = gpu_lerp(
            coil_sens[ix, iy1, iz, coil], coil_sens[ix1, iy1, iz, coil], wx,
        )
        c01 = gpu_lerp(
            coil_sens[ix, iy, iz1, coil], coil_sens[ix1, iy, iz1, coil], wx,
        )
        c11 = gpu_lerp(
            coil_sens[ix, iy1, iz1, coil], coil_sens[ix1, iy1, iz1, coil], wx,
        )
    end
    c0 = gpu_lerp(c00, c10, wy)
    c1 = gpu_lerp(c01, c11, wy)
    return gpu_lerp(c0, c1, wz)
end

# Scalar sensitivity methods called from Bloch GPU kernels.
@inline function KomaMRIBase.get_sens(receiver::BirdcageCoilGeometry, position, coil)
    x, y, z = position
    return KomaMRIBase.birdcage_wire_sensitivity(
        x, y, z, @inbounds(receiver.x[coil]), @inbounds(receiver.y[coil]), receiver.L,
    )
end

@inline function KomaMRIBase.get_sens(receiver::ArbitraryCoilSens, position, coil)
    x, y, z = position
    return gpu_interpolate_sensitivity(
        x, y, z, coil, receiver.x, receiver.y, receiver.z, receiver.coil_sens,
    )
end

# Whole-map materialization runs outside the kernel for static and custom receivers.
function KomaMRIBase.get_sens(receiver::ArbitraryCoilSens, x, y, z, ::KA.GPU)
    ncoils = get_n_coils(receiver)
    sens = similar(x, eltype(receiver.coil_sens), length(x), ncoils)
    coils = similar(x, Int, ncoils)
    coils .= axes(receiver.coil_sens, 4)
    sens .= gpu_interpolate_sensitivity.(
        reshape(x, :, 1),
        reshape(y, :, 1),
        reshape(z, :, 1),
        reshape(coils, 1, :),
        Ref(receiver.x),
        Ref(receiver.y),
        Ref(receiver.z),
        Ref(receiver.coil_sens),
    )
    return sens
end

# GPU simulation preallocation materializes static maps once. With phantom
# motion, the fixed birdcage/arbitrary models are evaluated at the moving spins.
prealloc_receiver(receiver::UniformCoilSens, _, _, ::NoMotion) = receiver
prealloc_receiver(
    receiver::UniformCoilSens, _, _, ::Union{Motion,MotionList},
) = receiver
prealloc_receiver(receiver, obj, backend, ::NoMotion) =
    get_sens(receiver, obj.x, obj.y, obj.z, backend)
prealloc_receiver(
    receiver::BirdcageCoilSens, obj, backend, ::Union{Motion,MotionList},
) = gpu_birdcage_geometry(receiver, eltype(obj.x), backend)
prealloc_receiver(
    receiver::ArbitraryCoilSens, _, _, ::Union{Motion,MotionList},
) = receiver

has_coil_sensitivities(::UniformCoilSens) = Val(false)
has_coil_sensitivities(_) = Val(true)
