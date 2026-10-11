module KomaMRICoreReactantExt

using Reactant
using KomaMRIBase: Phantom, DiscreteSequence, ArbitraryCoilSens, get_n_coils
import KomaMRIBase
using KomaMRICore
using ReactantCore: promote_to_traced
import KomaMRICore: _prepare_sequence,
    _run_simulation_block!, acquire_signal!, _sample_counter, gpu_grid_interval,
    gpu_interpolate_cell, gpu_coil_sample, prealloc_sensitivities, CoilSensitivities

function gpu_grid_interval(grid, query::Reactant.AnyTracedRArray)
    n = length(grid)
    low = Reactant.promote_to(Reactant.TracedRArray, ones(Int, size(query)))
    n == 1 && return low, zero.(query), query .== grid[1]
    high = low .* n
    for _ in 1:ceil(Int, log2(n))
        middle = (low .+ high) .÷ 2
        move_low = reshape(grid[vec(middle)], size(query)) .<= query
        low = ifelse.(move_low, middle, low)
        high = ifelse.(move_low, high, max.(middle, low .+ 1))
    end
    lower = reshape(grid[vec(low)], size(query))
    upper = reshape(grid[vec(high)], size(query))
    return low, (query .- lower) ./ (upper .- lower),
        (query .>= grid[1]) .& (query .<= grid[n])
end

function gpu_coil_sample(values::Reactant.AnyTracedRArray,
    ix::Reactant.AnyTracedRArray, iy, iz, coil)
    nx, ny, nz = size(values)[1:3]
    indices = @. ix + (iy - 1) * nx + (iz - 1) * nx * ny + (coil - 1) * nx * ny * nz
    return reshape(vec(values)[vec(indices)], size(indices))
end

function KomaMRIBase.get_sens(receiver::ArbitraryCoilSens,
    x::Reactant.AnyTracedRVector, y, z)
    traced(a) = Reactant.promote_to(Reactant.TracedRArray, a)
    x_grid, y_grid, z_grid = traced.((receiver.x, receiver.y, receiver.z))
    values = traced(receiver.coil_sens)
    coils = reshape(1:get_n_coils(receiver), 1, :)
    ix, wx, inside_x = gpu_grid_interval(x_grid, reshape(x, :, 1))
    iy, wy, inside_y = gpu_grid_interval(y_grid, reshape(y, :, 1))
    iz, wz, inside_z = gpu_grid_interval(z_grid, reshape(z, :, 1))
    values = gpu_interpolate_cell(ix, iy, iz, wx, wy, wz, coils,
        length(x_grid), length(y_grid), length(z_grid), values)
    return ifelse.(inside_x .& inside_y .& inside_z, values, zero(eltype(values)))
end

prealloc_sensitivities(receiver::ArbitraryCoilSens,
    obj::Phantom{<:Reactant.AnyTracedRVector}) =
    CoilSensitivities(KomaMRIBase.get_sens(receiver, obj.x, obj.y, obj.z), nothing)

function KomaMRIBase.get_sens!(sens::Reactant.AnyTracedRMatrix,
    receiver::ArbitraryCoilSens, x, y, z, _)
    sens .= KomaMRIBase.get_sens(receiver, x, y, z)
    return sens
end

_sample_counter(sig::Reactant.AnyTracedRArray) =
    promote_to_traced(firstindex(sig, 1))

function _prepare_sequence(::Phantom{<:Reactant.AnyTracedRVector}, seqd)
    traced(x) = Reactant.promote_to(Reactant.TracedRArray, x)
    DiscreteSequence(traced(seqd.Gx), traced(seqd.Gy), traced(seqd.Gz),
        traced(seqd.B1), traced(seqd.Δf), traced(seqd.ψ), seqd.ADC,
        seqd.excitation_bool, traced(seqd.t), traced(seqd.Δt))
end

function _run_simulation_block!(runner, obj::Phantom{<:Reactant.AnyTracedRVector},
    seqd, p, sig, samples, Xt, sys, method, groupsize, backend, preallocated;
    Nthreads)
    Nthreads == 1 || error("Compiled simulation requires Nthreads=1")
    block = seqd[p]
    signal = fill!(similar(sig, length(samples), size(sig)[2:end-1]...), 0)
    kernel = runner === KomaMRICore.run_spin_excitation_parallel! ?
        KomaMRICore.run_spin_excitation! : KomaMRICore.run_spin_precession!
    kernel(obj, block, signal, Xt, sys, method, groupsize, backend, preallocated)
    dims = ntuple(_ -> Colon(), ndims(sig) - 2)
    isempty(samples) || (sig[samples, dims..., 1] = signal)
    return nothing
end

function acquire_signal!(sig::Reactant.AnyTracedRMatrix, sample, Mxy, sensitivities, coords)
    signal = similar(sig, size(sig, 2))
    acquire_signal!(signal, Mxy, sensitivities, coords)
    sig[sample, :] = signal
    return nothing
end

end
