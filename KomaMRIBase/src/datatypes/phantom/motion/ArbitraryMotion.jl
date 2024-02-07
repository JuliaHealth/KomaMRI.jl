@with_kw mutable struct ArbitraryMotion{T<:Real} <: MotionModel
	etp_x
    etp_y
    etp_z
    etp_flags
end

function ArbitraryMotion( dur::AbstractVector{T},
                          K::Int,
                          Δx::AbstractArray{T, 2},
                          Δy::AbstractArray{T, 2},
                          Δz::AbstractArray{T, 2},
                          resetmag::BitMatrix ) where {T<:Real}

    etp(x) = extrapolate(x,Periodic())
    Ns = size(Δx)[1]
    limits = get_pieces_limits(dur,K)

    Δ = zeros(Ns,length(limits),3)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),Δx),1,length(dur)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),Δy),1,length(dur)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),Δz),1,length(dur)),zeros(Ns,1))

    itpx = sum(abs.(Δ[:,:,1]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,1], Gridded(Linear())) for i in 1:Ns] : nothing
    itpy = sum(abs.(Δ[:,:,2]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,2], Gridded(Linear())) for i in 1:Ns] : nothing
    itpz = sum(abs.(Δ[:,:,3]);dims=2) != zeros(Ns,1) ? [interpolate((limits,), Δ[i,:,3], Gridded(Linear())) for i in 1:Ns] : nothing
    flags =  any(resetmag) ? [interpolate((limits,), vcat(resetmag[i,:],Bool(0)), Gridded(Constant{Previous}())) for i in 1:Ns] : nothing

    etpx =      itpx  !== nothing ? map(etp, itpx)  : nothing 
    etpy =      itpy  !== nothing ? map(etp, itpy)  : nothing 
    etpz =      itpz  !== nothing ? map(etp, itpz)  : nothing  
    etpflags =  flags !== nothing ? map(etp, flags) : nothing

    ArbitraryMotion{Float64}(
        etp_x = etpx,
        etp_y = etpy,
        etp_z = etpz,
        etp_flags = etpflags
    )
end

ArbitraryMotion(etp_x, etp_y, etp_z, etp_flags) = ArbitraryMotion{eltype(eltype(@view(etp_x[1])))}(etp_x, etp_y, etp_z, etp_flags)



Base.getindex(motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}) = begin
    etp_x =     motion.etp_x !== nothing ?      motion.etp_x[p] :   nothing
    etp_y =     motion.etp_y !== nothing ?      motion.etp_y[p] :   nothing
    etp_z =     motion.etp_z !== nothing ?      motion.etp_z[p] :   nothing
    etp_flags = motion.etp_flags !== nothing ?  motion.etp_flags[p] : nothing
    ArbitraryMotion(
        etp_x = etp_x,
        etp_y = etp_y,
        etp_z = etp_z,
        etp_flags = etp_flags
    )
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
function get_pieces_limits(dur::AbstractVector, K::Int)
	steps = dur/K
	mat = reduce(hcat,[steps for i in 1:K])'
	limits = reshape(mat,(K*length(dur),))
	cumsum!(limits,limits)
	limits = vcat(0,limits)
    limits
end


function get_positions(motion::ArbitraryMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
    interpolate_spin_displacement(itp) = itp.(t)
    init = [x,y,z,0]
    positions = []
    
    for (i,field) in enumerate(fieldnames(ArbitraryMotion))
        etp = getproperty(motion,field)
        if etp !== nothing
            u = reduce(vcat,map(interpolate_spin_displacement,etp))
            xt = init[i] .+ u
            if size(xt,2) == 1
                xt = vec(xt)
            end
            push!(positions, xt)
        else
            if i == 4
                push!(positions,nothing)
            else
                push!(positions,init[i])
            end
        end
    end
    return positions
end


