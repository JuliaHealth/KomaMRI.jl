@with_kw mutable struct ArbitraryMotion{T} <: MotionModel
    # Segments
	dur::AbstractVector{T} = [1]
	K::Int = 2

	# Motion
	Δx::AbstractArray{T, 2} 
	Δy::AbstractArray{T, 2} = zeros(size(Δx))
	Δz::AbstractArray{T, 2} = zeros(size(Δx))

    resetmag::BitMatrix = zeros(size(Δx)[1],K) # (Ns x K)
end

function ArbitraryMotion(Ns::Int)
    ArbitraryMotion{Int}(Δx=zeros(Ns,1),
                         Δy=zeros(Ns,1),
                         Δz=zeros(Ns,1),
                         resetmag=zeros(Ns,2))       
end

export ArbitraryMotion

Base.getindex(mov::ArbitraryMotion, p::AbstractRange) = begin
    ArbitraryMotion(dur=mov.dur,
                    K=mov.K,
                    Δx=mov.Δx[p,:],
                    Δy=mov.Δy[p,:],
                    Δz=mov.Δz[p,:],
                    resetmag=mov.resetmag[p,:])
end

Base.getindex(mov::ArbitraryMotion, p::AbstractVector) = begin
    ArbitraryMotion(dur=mov.dur,
                    K=mov.K,
                    Δx=mov.Δx[p,:],
                    Δy=mov.Δy[p,:],
                    Δz=mov.Δz[p,:],
                    resetmag=mov.resetmag[p,:])
end

"""
    limits = get_pieces_limits(obj.mov)

Returns the pieces limits from dur and K values

Example: -----------------------------
    mov.dur = [1, 0.5]
    mov.K = 4

    limits = [0, 0.25, 0.5, 0.75, 1, 1.125, 1.25, 1.375, 1.5]
--------------------------------------
"""
function get_pieces_limits(mov::ArbitraryMotion)
	dur = mov.dur
	K   = mov.K

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
function get_itp_functions(mov::ArbitraryMotion{T}) where {T<:Real}
    Ns = size(mov.Δx)[1]
    limits = get_pieces_limits(mov)

    Δ = zeros(Ns,length(limits),3)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),mov.Δx),1,length(mov.dur)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),mov.Δy),1,length(mov.dur)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),mov.Δz),1,length(mov.dur)),zeros(Ns,1))

    itpx = sum(abs.(Δ[:,:,1]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,1], Gridded(Linear())) for i in 1:Ns] : nothing
    itpy = sum(abs.(Δ[:,:,2]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,2], Gridded(Linear())) for i in 1:Ns] : nothing
    itpz = sum(abs.(Δ[:,:,3]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,3], Gridded(Linear())) for i in 1:Ns] : nothing

    itp = [itpx, itpy, itpz]
    
    # dis = get_discontinuity(mov)
    flags = is_fluid(mov) ? [interpolate((limits,), vcat(mov.resetmag[i,:],Bool(0)), Gridded(Constant{Previous}())) for i in 1:Ns] : nothing

    itp,flags
end



"""
    Ux, Uy, Uz = initialize_motion(obj.mov, seqd.t)
"""
function initialize_motion(mov::ArbitraryMotion{T}, 
                           x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}; 
                           enable_gpu::Bool=false, gpu_device::Int=0, precision::AbstractString="f32") where {T<:Real}
    Ns = size(mov.Δx)[1]
    times = mod.(t,sum(mov.dur)) # Map time values between 0 and sum(dur)

    itp,flags = get_itp_functions(mov)
    is_mov_on = (itp .!== nothing) # tells in which dimensions we have movement
    is_fluid = (flags !== nothing)

    # Precision
    if precision == "f32"
        itp    = itp    |> f32
        flags  = flags  |> f32
    elseif precision == "f64"
        itp    = itp    |> f64
        flags  = flags  |> f64
    end
    # To GPU
    if enable_gpu
        device!(gpu_device)
        itp    = itp    |> gpu
        times  = times  |> gpu
        flags  = flags  |> gpu

        Ux = is_mov_on[1] ? reduce(hcat,[itp[1][i].(times) for i in 1:Ns])' : nothing
        Uy = is_mov_on[2] ? reduce(hcat,[itp[2][i].(times) for i in 1:Ns])' : nothing
        Uz = is_mov_on[3] ? reduce(hcat,[itp[3][i].(times) for i in 1:Ns])' : nothing

        resetmag = is_fluid ? Bool.(reduce(hcat,[flags[i].(times) for i in 1:Ns])') : nothing
    else
        Ux = is_mov_on[1] ? reduce(hcat,[itp[1][i](times) for i in 1:Ns])' : nothing
        Uy = is_mov_on[2] ? reduce(hcat,[itp[2][i](times) for i in 1:Ns])' : nothing
        Uz = is_mov_on[3] ? reduce(hcat,[itp[3][i](times) for i in 1:Ns])' : nothing

        resetmag = is_fluid ? Bool.(reduce(hcat,[flags[i](times) for i in 1:Ns])') : nothing
    end

    Ux, Uy, Uz, resetmag
end

function is_dynamic(mov::ArbitraryMotion{T}) where {T<:Real}
    itp = get_itp_functions(mov)
    return reduce(|,(itp .!== nothing))
end

function is_fluid(mov::ArbitraryMotion{T}) where {T<:Real}
    return any(mov.resetmag)
end

function get_discontinuity(mov::ArbitraryMotion)
    differences(x) = begin
        diff = zeros(size(x)[1],size(x)[2]+1)
        for j in 1:size(x)[1]
            diff[j,2:end-1] = abs.([x[j,i] - x[j,i-1] for i in range(2,size(x)[2])]) 
        end
        return diff
    end

    threshold = 0.5

    x = differences(mov.Δx) .> threshold
    y = differences(mov.Δy) .> threshold
    z = differences(mov.Δz) .> threshold

    hcat(repeat((x.|y.|z),1,length(mov.dur)),BitVector(zeros(size(x)[1])))
end
