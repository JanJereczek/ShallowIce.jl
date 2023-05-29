struct ComputationDomain{T<:AbstractFloat}
    L::T
    N::Int
    dx::T
    xH::Vector{T}
    xF::Vector{T}
    tspan::Tuple
    dt::Real
    t::Vector{T}
    bc::String
end

function ComputationDomain(L::T, N::Int, tspan::Tuple, dt::T,
    bc::String) where {T<:AbstractFloat}
    xH = range(0.0, stop = L, length = N)
    dx = xH[2] - xH[1]
    xF = xH .+ dx / 2
    t = collect(tspan[1]:dt:tspan[2])
    return ComputationDomain(L, N, dx, collect(xH), collect(xF), tspan, dt, t, bc)
end