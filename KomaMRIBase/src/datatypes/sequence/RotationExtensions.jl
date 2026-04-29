"""
    apply_rotations(seq; reverse=false)

Return a copy of `seq` with `QuaternionRot` extensions applied to the block
gradients. Rotation extensions are kept in `EXT`, so the operation can be
inverted and written back to Pulseq.
"""
function apply_rotations(seq; reverse=false)
    _has_rotation_extensions(seq.EXT) || return copy(seq)
    return _apply_rotations_to_owned_sequence(copy(seq); reverse)
end

function _apply_rotations_to_owned_sequence(seq; reverse=false)
    _has_rotation_extensions(seq.EXT) || return seq
    rotated_blocks = _rotated_blocks(seq.GR, seq.EXT; reverse)
    any(!isnothing, rotated_blocks) || return seq
    if _rotation_fits_gr_storage(seq.GR, rotated_blocks)
        _set_rotated_blocks!(seq.GR, rotated_blocks)
        return seq
    end
    GR = Matrix{Grad}(undef, size(seq.GR))
    GR .= seq.GR
    _set_rotated_blocks!(GR, rotated_blocks)
    return Sequence(GR, seq.RF, seq.ADC, seq.DUR, seq.EXT, seq.DEF)
end

_has_rotation_extensions(EXT) =
    any(block -> any(ext -> ext isa QuaternionRot, block), EXT)

function _rotated_block_gradients(gr, ext; reverse=false)
    out = gr
    has_rotation = false
    for e in ext
        e isa QuaternionRot || continue
        has_rotation = true
        out = (reverse ? e' : e) * out
    end
    return has_rotation ? out : nothing
end

function _rotated_blocks(GR, EXT; reverse=false)
    rotated_blocks = Vector{Union{Nothing,AbstractVector{<:Grad}}}(nothing, size(GR, 2))
    Threads.@threads for b in eachindex(EXT)
        rotated_blocks[b] = _rotated_block_gradients(view(GR, :, b), EXT[b]; reverse)
    end
    return rotated_blocks
end

function _apply_rotation_extensions_to_gradients!(GR, EXT; reverse=false)
    rotated_blocks = _rotated_blocks(GR, EXT; reverse)
    any(!isnothing, rotated_blocks) || return GR
    _rotation_fits_gr_storage(GR, rotated_blocks) && return _set_rotated_blocks!(GR, rotated_blocks)

    out = Matrix{Grad}(undef, size(GR))
    out .= GR
    _set_rotated_blocks!(out, rotated_blocks)
    return _concrete_event_array_copy(out)
end

function _rotation_fits_gr_storage(GR, rotated_blocks)
    T = eltype(GR)
    for gr in rotated_blocks
        isnothing(gr) && continue
        all(g -> g isa T, gr) || return false
    end
    return true
end

function _set_rotated_blocks!(GR, rotated_blocks)
    Threads.@threads for b in eachindex(rotated_blocks)
        gr = rotated_blocks[b]
        isnothing(gr) && continue
        view(GR, :, b) .= gr
    end
    return GR
end
