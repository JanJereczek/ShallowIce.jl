mutable struct IcesheetState{T<:AbstractFloat}
    h::Vector{T}        # thickness
    b::Vector{T}        # bedrock
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
    return IcesheetState(h, b, alpha, beta, gamma, delta, f, g)
end

struct SuperStruct
    p::Params
    omega::ComputationDomain
    iss::IcesheetState
end

"""
    right_mean(x::Vector, N::Int)
Return mean on staggered grid, given by:

"""
function right_mean(x::Vector, N::Int)
    return (view(x, 1:N-1) + view(x, 2:N)) ./ 2
end

function forward_sia(sstruct::SuperStruct)
    nt = Int(sstruct.p.tspan[2] ÷ sstruct.omega.dt)
    ht = zeros(sstruct.p.N, nt)
    for k in 1:nt
        forwardstep_sia!(sstruct, sstruct.omega.dt)
        ht[:, k] .= copy(sstruct.iss.h)
    end
    return ht
end

function forwardstep_sia!(
    sstruct::SuperStruct,
    dt::T,
) where {T<:AbstractFloat}

    omega, p, iss = sstruct.omega, sstruct.p, sstruct.iss
    N, dx = omega.N, omega.dx
    n = p.n
    h, b = iss.h, iss.b
    dtdx = dt / omega.dx ^ 2

    alpha, beta, gamma, delta = iss.alpha, iss.beta, iss.gamma, iss.delta
    f, g = iss.f, iss.g

    dsdx = fdx( h + b, dx, N )
    d = sstruct.p.fd .* abs.(dsdx) .^ (n - 1) .* right_mean(h, N) .^ (n + 2)
    a = compute_mass_balance(N)

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
    return nothing
end

function compute_mass_balance(N::Int)
    return fill(0.3, N)
end

function gaussian_mean(x::Vector)


function normed_gaussian(x, mu, sigma)
    gaussian = exp( (x-mu)^2 / sigma )
    return gaussian ./ maximum(gaussian)
end

function dbdt!(dbdt, b, sstruct, t)
    return 1/p.tau .* (sstruct.iss.b .-
        (sstruct.rho_ice / sstruct.rho_mantle) .*
        gaussian_mean(sstruct.iss.h))
end