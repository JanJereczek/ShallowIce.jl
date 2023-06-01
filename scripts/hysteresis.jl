include("intro.jl")
using JLD2

function a_shift(
    t::Real,
    amax::Real = 0.0025;
    t_a_up::Real = 900e3,
    t_a_const::Real = 1000e3,
    t_a_down::Real = 1900e3,
)
    m_up = amax / t_a_up
    m_down = -m_up
    if t < t_a_up
        return m_up * t
    elseif t < t_a_const
        return amax
    elseif t < t_a_down
        return m_down * (t-t_a_down)
    else
        return 0.0
    end
end

# t = 0:100:2e6
# lines(a_shift.(t))

function main(; N = 101, tend = 2e6, dt = 2.0, dt_out = 1000.0)
    L = 1.5e6
    tspan = (0.0, tend)
    bc = "zero_flow"

    omega = ComputationDomain(L, N, tspan, dt, bc)
    bedhead = 2e3               # m
    bedslope = -0.002           # (m/m)
    bedbump_center = 1e6        # m
    bedbump_height = 500.0      # m
    bedbump_width = 100e3       # m
    h0 = fill(0.0, N)
    b0 = bumpy_bed(omega.xH, bedhead, bedslope, bedbump_center, bedbump_height, bedbump_width)

    a1 = 0.25
    a2 = -0.4
    init_a1, init_a2 = 0.0, a2-a1
    shifting_linear_accumulation(t, x) = a_shift(t, a1) .+
        linear_accumulation(t, x, a1 = init_a1, a2 = init_a2)

    p = Params(isostasy_on = true, accumulation = shifting_linear_accumulation)
    iss = IcesheetState(N, h0, b0)
    sstruct = SuperStruct(p, omega, iss)
    ht, bt = forward_sia(sstruct, dt_out = dt_out)
    t_out = Int.(collect(tspan[1]:dt_out:tspan[2]))
    a_out = [shifting_linear_accumulation(t, omega.xH) for t in t_out]
    jldsave("data/hysteresis_N=$(N)_tend=$(tend).jld2"; ht, bt, t_out, a_out)
end

main(N = 51, dt = 10.0)