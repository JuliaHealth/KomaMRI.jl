@with_kw mutable struct SimpleMotion <: MotionModel
    ux::Function = (x,y,z,t)->0
	uy::Function = (x,y,z,t)->0
	uz::Function = (x,y,z,t)->0
end
export SimpleMotion

Base.getindex(mov::SimpleMotion, p::AbstractRange) = mov
Base.getindex(mov::SimpleMotion, p::AbstractVector) = mov

"""
    Ux, Uy, Uz = initialize_motion(obj.mov, seqd.t)
"""
function initialize_motion(mov::SimpleMotion, 
                           x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}; 
                           enable_gpu::Bool=false, gpu_device::Int=0, precision::AbstractString="f32") where {T<:Real}
 
    Ns = length(x)
                           
    Ux = zeros(Ns,1) .+ mov.ux(x, y, z, t')
    Uy = zeros(Ns,1) .+ mov.uy(x, y, z, t')
    Uz = zeros(Ns,1) .+ mov.uz(x, y, z, t')

    Ux = sum(abs.(Ux);dims=2) != zeros(Ns,1) ? Ux : nothing
    Uy = sum(abs.(Uy);dims=2) != zeros(Ns,1) ? Uy : nothing
    Uz = sum(abs.(Uz);dims=2) != zeros(Ns,1) ? Uz : nothing

    # Precision
    if precision == "f32"
        Ux  = Ux  |> f32
        Uy  = Uy  |> f32
        Uz  = Uz  |> f32
    elseif precision == "f64"
        Ux  = Ux  |> f64
        Uy  = Uy  |> f64
        Uz  = Uz  |> f64
    end
    # To GPU
    if enable_gpu
        device!(gpu_device)
        Ux  = Ux  |> gpu
        Uy  = Uy  |> gpu
        Uz  = Uz  |> gpu
    end

    Ux, Uy, Uz
end