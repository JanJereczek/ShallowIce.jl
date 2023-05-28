struct ComputationDomain{T<:AbstractFloat}
    L::T
    N::Int
    dx::T
    xH::Vector{T}
    xF::Vector{T}
    tspan::Tuple
end

function ComputationDomain(L::T, N::Int, tspan::Tuple) where {T<:AbstractFloat}
    xH = range(0.0, stop = L, length = N)
    dx = xH[2] - xH[1]
    xF = xH .+ dx / 2
    return ComputationDomain(L, N, dx, collect(xH), collect(xF), tspan)
end