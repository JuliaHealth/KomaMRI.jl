@with_kw mutable struct ArbitraryMotion{T} <: MotionModel{T}
	dur::AbstractVector{T} = [1.0]
	K::Int = 2

	# Motion
	Δx::AbstractArray{T, 2} 
	Δy::AbstractArray{T, 2} = zeros(size(Δx))
	Δz::AbstractArray{T, 2} = zeros(size(Δx))

    resetmag::BitMatrix = zeros(size(Δx)[1],K) # (Ns x K)
end

function ArbitraryMotion(Ns::Int)
    ArbitraryMotion(Δx=zeros(Ns,1))       
end

include("ExplicitArbitraryMotion.jl")


Base.getindex(motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}) = begin
    ArbitraryMotion(dur=motion.dur,
                    K=motion.K,
                    Δx=motion.Δx[p,:],
                    Δy=motion.Δy[p,:],
                    Δz=motion.Δz[p,:],
                    resetmag=motion.resetmag[p,:])
end
Base.getindex(motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                       q::Union{AbstractRange,AbstractVector,Colon}) = motion[p]


+(m1::ArbitraryMotion,m2::ArbitraryMotion) = begin
    if m1.K == m2.K
        return ArbitraryMotion(dur=m1.dur,
                               K = m1.K,
                               Δx = vcat(m1.Δx,m2.Δx),
                               Δy = vcat(m1.Δy,m2.Δy),
                               Δz = vcat(m1.Δz,m2.Δz),
                               resetmag = vcat(m1.resetmag,m2.resetmag))
    else
        return SimpleMotion()
    end
end

"""
    limits = get_pieces_limits(obj.motion)

Returns the pieces limits from dur and K values

Example: -----------------------------
    motion.dur = [1, 0.5]
    motion.K = 4

    limits = [0, 0.25, 0.5, 0.75, 1, 1.125, 1.25, 1.375, 1.5]
--------------------------------------
"""
function get_pieces_limits(motion::ArbitraryMotion)
	dur = motion.dur
	K   = motion.K

	steps = dur/K
	mat = reduce(hcat,[steps for i in 1:K])'
	limits = reshape(mat,(K*length(dur),))
	cumsum!(limits,limits)
	limits = vcat(0,limits)
    limits
end



"""
    itp = get_itp_functions(obj)

Returns an array of motion interpolation functions from a phantom
"""
function get_itp_functions(motion::ArbitraryMotion)
    Ns = size(motion.Δx)[1]
    limits = get_pieces_limits(motion)

    Δ = zeros(Ns,length(limits),3)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),motion.Δx),1,length(motion.dur)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),motion.Δy),1,length(motion.dur)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),motion.Δz),1,length(motion.dur)),zeros(Ns,1))

    itpx = sum(abs.(Δ[:,:,1]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,1], Gridded(Linear())) for i in 1:Ns] : nothing
    itpy = sum(abs.(Δ[:,:,2]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,2], Gridded(Linear())) for i in 1:Ns] : nothing
    itpz = sum(abs.(Δ[:,:,3]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,3], Gridded(Linear())) for i in 1:Ns] : nothing
    

    flags = is_fluid(motion) ? [interpolate((limits,), vcat(motion.resetmag[i,:],Bool(0)), Gridded(Constant{Previous}())) for i in 1:Ns] : nothing

    [itpx, itpy, itpz, flags]
end



"""
    Ux, Uy, Uz = initialize_motion(obj.motion, seqd.t)
"""
function initialize_motion(motion::ArbitraryMotion, t::AbstractVector{T}, sim_params::Dict) where {T<:Real}
    enable_gpu = sim_params["enable_gpu"]
    gpu_device = sim_params["gpu_device"]
    precision  = sim_params["precision"]

    Ns = size(motion.Δx)[1]
    times = mod.(t,sum(motion.dur)) # Map time values between 0 and sum(dur)

    itp = get_itp_functions(motion)

    U = [itp[i] !== nothing ? zeros(Ns,length(times)) : nothing for i in 1:4]

    if reduce(|,itp .!== nothing)
        # Precision
        if precision == "f32"
            itp    = itp    |> f32
            U      = U      |> f32
        elseif precision == "f64"
            itp    = itp    |> f64
            U      = U      |> f64
        end

        # To GPU
        if enable_gpu
            device!(gpu_device)

            times  = times  |> gpu

            part_size = 1e4
            Nparts = Int(ceil(Ns/part_size))
            parts = kfoldperm(Ns, Nparts; type="ordered")

            progress_bar = Progress(Nparts*length(U);desc="Initializing motion...")
        
            for (i, u) = enumerate(U)
                for (part, p) = enumerate(parts)
                    interp = itp[i] !== nothing ? itp[i][p]  |> gpu : nothing
                    mat    =      u !== nothing ? u[p,:]     |> gpu : nothing

                    if u !== nothing
                        interpolate_mov!(mat, interp, times)
                        U[i][p,:] = mat |> cpu
                    end
                    next!(progress_bar)
                end
            end 

        # In CPU
        else 
            progress_bar = Progress(length(U);desc="Initializing motion...")
            for (i, u) = enumerate(U)
                u = itp[i] !== nothing ? reduce(hcat,[itp[i][j](times) for j in 1:Ns])' : nothing
                #Update progress
                next!(progress_bar)
            end
        end
        return ExplicitArbitraryMotion(Ux    = U[1],
                                       Uy    = U[2],
                                       Uz    = U[3],
                                       flags = U[4])
    else
        return NoMotion()
    end
end


function interpolate_mov!(u, itp, times)
    if itp !== nothing
        u .= reduce(hcat,[func.(times) for func in itp])'
    end
end


function is_dynamic(motion::ArbitraryMotion)
    itp = get_itp_functions(motion)
    return reduce(|,(itp[1:3] .!== nothing))
end

function is_fluid(motion::ArbitraryMotion) 
    return any(motion.resetmag)
end


# UNUSED
function get_discontinuity(motion::ArbitraryMotion)
    differences(x) = begin
        diff = zeros(size(x)[1],size(x)[2]+1)
        for j in 1:size(x)[1]
            diff[j,2:end-1] = abs.([x[j,i] - x[j,i-1] for i in range(2,size(x)[2])]) 
        end
        return diff
    end

    threshold = 0.5

    x = differences(motion.Δx) .> threshold
    y = differences(motion.Δy) .> threshold
    z = differences(motion.Δz) .> threshold

    hcat(repeat((x.|y.|z),1,length(motion.dur)),BitVector(zeros(size(x)[1])))
end



