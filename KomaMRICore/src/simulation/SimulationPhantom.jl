struct SimulationPhantom{
    T<:Real,
    X<:AbstractVector{T},
    Y<:AbstractVector{T},
    Z<:AbstractVector{T},
    R<:AbstractVector{T},
    T1A<:AbstractVector{T},
    T2A<:AbstractVector{T},
    WA<:AbstractVector{T},
    M,
} <: AbstractPhantom{T}
    x::X
    y::Y
    z::Z
    ρ::R
    T1::T1A
    T2::T2A
    Δw::WA
    motion::M
end

@functor SimulationPhantom

@inline KomaMRIBase.get_x(obj::SimulationPhantom) = obj.x
@inline KomaMRIBase.get_y(obj::SimulationPhantom) = obj.y
@inline KomaMRIBase.get_z(obj::SimulationPhantom) = obj.z
@inline KomaMRIBase.get_ρ(obj::SimulationPhantom) = obj.ρ
@inline KomaMRIBase.get_T1(obj::SimulationPhantom) = obj.T1
@inline KomaMRIBase.get_T2(obj::SimulationPhantom) = obj.T2
@inline KomaMRIBase.get_Δw(obj::SimulationPhantom) = obj.Δw
@inline KomaMRIBase.get_motion(obj::SimulationPhantom) = obj.motion

Base.length(obj::SimulationPhantom) = length(obj.ρ)
Base.size(obj::SimulationPhantom) = size(obj.ρ)

function Base.view(obj::SimulationPhantom, p)
    return SimulationPhantom(
        @view(obj.x[p]),
        @view(obj.y[p]),
        @view(obj.z[p]),
        @view(obj.ρ[p]),
        @view(obj.T1[p]),
        @view(obj.T2[p]),
        @view(obj.Δw[p]),
        @view(obj.motion[p]),
    )
end

_convert_property(::Type{T}, x::AbstractVector{T}) where {T<:Real} = x
_convert_property(::Type{T}, x::AbstractVector) where {T<:Real} = T.(x)

function SimulationPhantom(obj::AbstractPhantom, ::Type{T}) where {T<:Real}
    x = get_x(obj)
    y = get_y(obj)
    z = get_z(obj)
    ρ = get_ρ(obj)
    T1 = get_T1(obj)
    T2 = get_T2(obj)
    Δw = get_Δw(obj)
    N = length(ρ)
    for (name, property) in
        ((:x, x), (:y, y), (:z, z), (:T1, T1), (:T2, T2), (:Δw, Δw))
        length(property) == N ||
            throw(DimensionMismatch("Phantom property `$name` has length $(length(property)); expected $N."))
    end
    return SimulationPhantom(
        _convert_property(T, x),
        _convert_property(T, y),
        _convert_property(T, z),
        _convert_property(T, ρ),
        _convert_property(T, T1),
        _convert_property(T, T2),
        _convert_property(T, Δw),
        paramtype(T, get_motion(obj)),
    )
end
