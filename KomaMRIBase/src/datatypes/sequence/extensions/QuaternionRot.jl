mutable struct QuaternionRot <: Extension
    q0::Float64
    qx::Float64
    qy::Float64
    qz::Float64
end

get_scale(::Type{QuaternionRot}) = [1.0 1.0 1.0 1.0]
get_scanf_format(::Type{QuaternionRot}) = "%f %f %f %f"
get_EXT_type_from_symbol(::Val{:ROTATIONS}) = QuaternionRot
Base.transpose(q::QuaternionRot) = QuaternionRot(q.q0, -q.qx, -q.qy, -q.qz)
Base.adjoint(q::QuaternionRot) = transpose(q)

function Base.:*(q::QuaternionRot, grads)
    @assert length(grads) == 3 "A sequence block must contain x, y and z gradient channels."

    qnorm = sqrt(abs2(q.q0) + abs2(q.qx) + abs2(q.qy) + abs2(q.qz))
    qnorm > 0 || throw(ArgumentError("Rotation quaternion must have non-zero norm."))

    q0, qx, qy, qz = q.q0 / qnorm, q.qx / qnorm, q.qy / qnorm, q.qz / qnorm

    rot = [
        (1 - 2 * (qy^2 + qz^2))   2 * (qx * qy - q0 * qz)   2 * (qx * qz + q0 * qy)
        2 * (qx * qy + q0 * qz)   (1 - 2 * (qx^2 + qz^2))   2 * (qy * qz - q0 * qx)
        2 * (qx * qz - q0 * qy)   2 * (qy * qz + q0 * qx)   (1 - 2 * (qx^2 + qy^2))
    ]
    return rot * grads
end

function apply_rotations!(seq; reverse=false)
    for b in eachindex(seq.EXT)
        block_ext = seq.EXT[b]
        isempty(block_ext) && continue

        for ext in block_ext
            ext isa QuaternionRot || continue
            seq.GR[:, b] = (reverse ? ext' : ext) * seq.GR[:, b]
        end
    end
    return seq
end
