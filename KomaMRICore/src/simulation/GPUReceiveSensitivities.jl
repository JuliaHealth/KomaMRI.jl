receiver_sensitivities(receiver, x, y, z, ::KA.CPU) = get_sens(receiver, x, y, z)
receiver_sensitivities(receiver, x, y, z, ::KA.GPU) = get_sens(receiver, x, y, z)
receiver_sensitivities(receiver, x, y, z) =
    receiver_sensitivities(receiver, x, y, z, KA.get_backend(x))

@inline function gpu_grid_interval(grid, query)
    n = length(grid)
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
    c0 = gpu_lerp(c00, c10, wy)
    c1 = gpu_lerp(c01, c11, wy)
    return gpu_lerp(c0, c1, wz)
end

function receiver_sensitivities(receiver::ArbitraryCoilSens, x, y, z, ::KA.GPU)
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
