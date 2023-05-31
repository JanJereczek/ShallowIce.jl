mutable struct IcesheetState{T<:AbstractFloat}
    h::Vector{T}        # thickness
    b::Vector{T}        # bedrock
    h0::Vector{T}
    b0::Vector{T}
    alpha::Vector{T}    # for semi-implicit scheme
    beta::Vector{T}     # for semi-implicit scheme
    gamma::Vector{T}    # for semi-implicit scheme
    delta::Vector{T}    # for semi-implicit scheme
    f::Vector{T}        # for semi-implicit scheme
    g::Vector{T}        # for semi-implicit scheme
end

function IcesheetState(
    N::Int,
    h = fill(0.0, N),
    b = fill(0.0, N);
    alpha = fill(0.0, N),
    beta = fill(0.0, N),
    gamma = fill(0.0, N),
    delta = fill(0.0, N),
    f = fill(0.0, N),
    g = fill(0.0, N),
)
    return IcesheetState(h, b, copy(h), copy(b), alpha, beta, gamma, delta, f, g)
end

struct SuperStruct
    p::Params
    omega::ComputationDomain
    iss::IcesheetState
end

"""
    right_mean(x::Vector, N::Int)

Return (right-sided) mean on staggered grid.
"""
function right_mean(x::Vector, N::Int)
    return (view(x, 1:N-1) + view(x, 2:N)) ./ 2
end

function forward_sia(sstruct::SuperStruct) # ; saving_stride::Int = 10
    nt = Int(sstruct.omega.tspan[2] ÷ sstruct.omega.dt)
    ht = zeros(sstruct.omega.N, nt)
    bt = similar(ht)
    for k in 1:nt
        forwardstep_sia!(sstruct, sstruct.omega.dt)
        ht[:, k] .= copy(sstruct.iss.h)
        bt[:, k] .= copy(sstruct.iss.b)
    end
    return ht, bt
end

function forwardstep_sia!(
    sstruct::SuperStruct,
    dt::T,
) where {T<:AbstractFloat}

    omega, p, iss = sstruct.omega, sstruct.p, sstruct.iss
    N, dx = omega.N, omega.dx
    n = p.n
    dtdx = dt / omega.dx ^ 2
    h, b = iss.h, iss.b

    alpha, beta, gamma, delta = iss.alpha, iss.beta, iss.gamma, iss.delta
    f, g = iss.f, iss.g

    dsdx = fdx( h + b, dx, N )
    a_raw = p.accumulation(omega.xH)
    a = a_raw .* ((h .> 0) .| (a_raw .> 0))

    d = sstruct.p.fd .* abs.(dsdx) .^ (n - 1) .* right_mean(h, N) .^ (n + 2)
    if omega.bc == "zero_flow"
        d[1] = 0.0
        d[end] = 0.0
    end

    if omega.bc == "zero_flow"
        sstruct.iss.alpha[1] = d[2] * dtdx
        sstruct.iss.gamma[1] = d[1] * dtdx
        sstruct.iss.beta[1] = 1 .+ gamma[1] .+ alpha[1]
        sstruct.iss.delta[1] = h[1] + a[1]*dt + alpha[1] * b[2] -
            (beta[1]-1)*b[1] + gamma[1]*b[2]
        sstruct.iss.f[1] = gamma[1] / (beta[1] - alpha[1] * f[2])
        sstruct.iss.g[1] = (delta[1] + alpha[1] * g[2]) /
            (beta[1] - alpha[1] * f[2])
        sstruct.iss.h[1] = g[1] + f[1] * h[2]
    end

    for j in 2:N-1
        sstruct.iss.alpha[j] = d[j-1] * dtdx
        sstruct.iss.gamma[j] = d[j] * dtdx
        sstruct.iss.beta[j] = 1 .+ gamma[j] .+ alpha[j]
        sstruct.iss.delta[j] = h[j] + a[j]*dt + alpha[j] * b[j-1] -
            (beta[j]-1)*b[j] + gamma[j]*b[j+1]
        sstruct.iss.f[j] = gamma[j] / (beta[j] - alpha[j] * f[j-1])
        sstruct.iss.g[j] = (delta[j] + alpha[j] * g[j-1]) /
            (beta[j] - alpha[j] * f[j-1])
    end
    for j in N-1:-1:1
        sstruct.iss.h[j] = g[j] + f[j] * h[j+1]
    end

    if p.isostasy_on
        forwardstep_isostasy!(sstruct)
    end

    return nothing
end

function forwardstep_isostasy!(sstruct)
    omega, p, iss = sstruct.omega, sstruct.p, sstruct.iss
    sstruct.iss.b .+= omega.dt/p.tau .* (iss.b0 .- iss.b .- (p.rho_ice / p.rho_mantle) .*
        gaussian_filter(omega.xH, iss.h) )
    return nothing
end