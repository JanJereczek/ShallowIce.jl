mutable struct IcesheetState{T<:AbstractFloat}
    t::T
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
    return IcesheetState(0.0, copy(h), copy(b), copy(h), copy(b),
        alpha, beta, gamma, delta, f, g)
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

function forward_sia(sstruct::SuperStruct; dt_out::Real = 10.0)
    nt = Int(sstruct.omega.tspan[2] ÷ sstruct.omega.dt)
    nt_out = Int(sstruct.omega.tspan[2] ÷ dt_out)

    ht_out = zeros(sstruct.omega.N, nt_out)
    bt_out = similar(ht_out)
    k_out = 1

    for k in 1:nt
        forwardstep_sia!(sstruct)
        if k*sstruct.omega.dt > k_out*dt_out
            ht_out[:, k_out] .= copy(sstruct.iss.h)
            bt_out[:, k_out] .= copy(sstruct.iss.b)
            println("t = $(Int(k_out*dt_out))")
            k_out += 1
        end
    end
    return ht_out, bt_out
end

function forwardstep_sia!(sstruct::SuperStruct)

    # Extract fields from sstruct
    omega, p, iss = sstruct.omega, sstruct.p, sstruct.iss
    N, dx, dt = omega.N, omega.dx, omega.dt
    dtdx = dt / omega.dx ^ 2
    n = p.n
    h, b = iss.h, iss.b
    alpha, beta, gamma, delta = iss.alpha, iss.beta, iss.gamma, iss.delta
    f, g = iss.f, iss.g

    # Compute surface slope and mass balance
    dsdx = fdx( h + b, dx, N )
    a_raw = p.accumulation(iss.t, omega.xH)
    a = a_raw .* ((h .> 0) .| (a_raw .> 0))

    # Compute diffusion
    d = sstruct.p.fd .* abs.(dsdx) .^ (n - 1) .* right_mean(h, N) .^ (n + 2)
    if omega.bc == "zero_flow"
        d[1] = 0.0
        d[end] = 0.0
    end

    # If zero flow, assume that cell j=0 is same as j=2
    # (and divide is at j = 1, first point of the domain).
    # Based on this, compute the semi-implicit integration.
    if omega.bc == "zero_flow"
        sstruct.iss.alpha[1] = d[2] * dtdx
        sstruct.iss.gamma[1] = d[1] * dtdx
        sstruct.iss.beta[1] = 1 .+ gamma[1] .+ alpha[1]
        sstruct.iss.delta[1] = h[1] + a[1]*dt + alpha[1] * b[2] -
            (beta[1]-1)*b[1] + gamma[1]*b[2]
        sstruct.iss.f[1] = gamma[1] / (beta[1] - alpha[1] * f[2])
        sstruct.iss.g[1] = (delta[1] + alpha[1] * g[2]) /
            (beta[1] - alpha[1] * f[2])
        # sstruct.iss.h[1] = g[1] + f[1] * h[2]
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

    # Update isostaic adjustment and go to next time step.
    if p.isostasy_on
        forwardstep_isostasy!(sstruct)
    end
    sstruct.iss.t += dt

    return nothing
end

function forwardstep_isostasy!(sstruct)
    omega, p, iss = sstruct.omega, sstruct.p, sstruct.iss
    sstruct.iss.b .+= omega.dt/p.tau .* (iss.b0 .- iss.b .- (p.rho_ice / p.rho_mantle) .*
        gaussian_filter(omega.xH, iss.h) )
    return nothing
end

# function implicit_integration!(sstruct, d, jm1, j, jp1, dt, dtdx)
#     sstruct.iss.alpha[j] = d[jm1] * dtdx
#     sstruct.iss.gamma[j] = d[j] * dtdx
#     sstruct.iss.beta[j] = 1 .+ sstruct.iss.gamma[j] .+ sstruct.iss.alpha[j]
#     sstruct.iss.delta[j] = h[j] + a[j]*dt + alpha[j] * b[j-1] -
#         (beta[j]-1)*b[j] + gamma[j]*b[j+1]
#     sstruct.iss.f[j] = gamma[j] / (beta[j] - alpha[j] * f[j-1])
#     sstruct.iss.g[j] = (delta[j] + alpha[j] * g[j-1]) /
#         (beta[j] - alpha[j] * f[j-1])
# end