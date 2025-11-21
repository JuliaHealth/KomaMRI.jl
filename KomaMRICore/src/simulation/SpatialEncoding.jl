

@inline function get_Bz(seq::DiscreteSequence, x, y, z, s_idx)
    return x .* seq.Gx[s_idx] .+ y .* seq.Gy[s_idx] .+ z .* seq.Gz[s_idx]
end
