include("intro.jl")
using JLD2

function linear_smb(z::Vector; zeq = 1000, c = 0.025 / 3000)
    return c .* (z.-zeq)
end

function zeq_shift(
    t::Real;
    zeq_min::Real = 1900,
    zeq_max::Real = 2100,
    t_zeq_down::Real = 900e3,
    t_zeq_const::Real = 1000e3,
    t_zeq_up::Real = 1900e3,
)
    dz = zeq_max - zeq_min
    m = dz / t_zeq_down
    if t < t_zeq_down
        return zeq_max - m * t
    elseif t < t_zeq_const
        return zeq_min
    elseif t < t_zeq_up
        return zeq_min .+ m * (t-t_zeq_const)
    else
        return zeq_max
    end
end
t = 0:100:2e6
using CairoMakie
lines(zeq_shift.(t))
hyst_smb(sstruct::SuperStruct) = hyst_smb(sstruct.iss.t, sstruct.iss.b + sstruct.iss.h)
hyst_smb(t, z) = linear_smb(z, zeq = zeq_shift(t))

function main(; N = 101, tend = 2e6, dt = 2.0, dt_out = 1000.0, isostasy_on = false)
    L = 1.5e6  
    tspan = (0.0, tend)
    bc = "zero_flow"
    t_out = Int.(collect(tspan[1]:dt_out:tspan[2]))

    omega = ComputationDomain(L, N, tspan, dt, bc)
    bedhead = 2e3               # m
    bedslope = -0.002           # (m/m)
    bedbump_center = 1e6        # m
    bedbump_height = 500.0      # m
    bedbump_width = 100e3       # m
    h0 = fill(0.0, N)
    b0 = bumpy_bed(omega.xH, bedhead, bedslope, bedbump_center, bedbump_height, bedbump_width)

    p = Params(accumulation = hyst_smb, isostasy_on = isostasy_on)
    iss = IcesheetState(N, h0, b0)
    sstruct = SuperStruct(p, omega, iss)
    ht, bt = forward_sia(sstruct, dt_out = dt_out)
    t_out = Int.(collect(tspan[1]:dt_out:tspan[2]))
    z_out = zeq_shift.(t_out)
    jldsave("data/hysteresis_N=$(N)_tend=$(tend)_isostasy=$(isostasy_on).jld2";
        ht, bt, t_out, z_out)
end

main(N = 76, dt = 5.0)