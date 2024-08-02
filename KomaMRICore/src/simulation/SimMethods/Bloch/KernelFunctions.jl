using KernelAbstractions: @kernel, @Const, @index, @uniform, @groupsize, @localmem

#=
macro myunroll(N)
    exp = []
    var1 = esc(quote
        for t=1:$N
            β[i_g] = -im * nxy[i_g, t] * sin(φ[i_g, t] / 2)
        end
    end)
    return var1
    push!(exp, var1)
    return esc(quote
        $exp[1]
        $exp[2]
    end)
end
=#

@kernel function apply_excitation!(α, β, Mxy_tmp, Mz_tmp, Mxy, Mz, Bz, B, φ, nxy, nz, ΔT1, ΔT2, ρ)
    i_g = @index(Global)
    #i_l = @index(Local)

    T = eltype(Bz)
    #@uniform N = @groupsize()[1]
    #@uniform N_Δt = size(φ, 2)

    #s_α = @localmem Complex{T} (N,)
    #s_β = @localmem Complex{T} (N,)
    #s_Mxy_tmp = @localmem Complex{T} (N,)
    #s_Mz_tmp = @localmem T (N,)

    N = size(φ, 2)
    #KA.Extras.LoopInfo.@unroll 
    @inbounds for t=1:N
        α[i_g] = cos(φ[i_g, t] / 2) - 1im * nz[i_g, t] * sin(φ[i_g, t] / 2)
        β[i_g] = -im * nxy[i_g, t] * sin(φ[i_g, t] / 2)
        Mxy_tmp[i_g] = 2 * conj(α[i_g]) * β[i_g] * Mz[i_g] + conj(α[i_g])^2 * Mxy[i_g] - β[i_g]^2 * conj(Mxy[i_g])
        Mz_tmp[i_g] = (abs(α[i_g])^2 - abs(β[i_g])^2) * Mz[i_g] - 2 * real(α[i_g] * β[i_g] * conj(Mxy[i_g]))
        Mxy[i_g] = Mxy_tmp[i_g] * exp(-ΔT2[i_g, t])
        Mz[i_g] = Mz_tmp[i_g] * exp(-ΔT1[i_g, t]) + ρ[i_g] * (1 - exp(-ΔT1[i_g, t]))
    end
    #@myunroll n
end
