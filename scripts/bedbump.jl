include("intro.jl")
using CairoMakie

N = 101
L = 1.5e6
tspan = (0.0, 20e3)
dt = 1.0
bc = "zero_flow"

omega = ComputationDomain(L, N, tspan, dt, bc)
bedhead = 3e3               # m
bedslope = -0.003           # (m/m)
bedbump_center = 1e6        # m
bedbump_height = 500.0      # m
bedbump_width = 100e3       # m

h0 = fill(0.0, N)
b0 = bumpy_bed(omega.xH, bedhead, bedslope, bedbump_center, bedbump_height, bedbump_width)
lines(b0)

p = Params(isostasy_on = true, accumulation = linear_accumulation)
iss = IcesheetState(N, h0, b0)

sstruct = SuperStruct(p, omega, iss)
ht, bt = forward_sia(sstruct)

idx = Observable(1)
stride = 100
hplot = @lift(ht[:, $idx])
bplot = @lift(bt[:, $idx])
splot = @lift((bt .+ ht)[:, $idx])
fig, ax, lh = lines(splot)
lb = lines!(ax, bplot)
ylims!(ax, (-1500, 4000))
record(fig, "plots/bumpy_isostasy_N=$(N)_$bc.gif", axes(ht, 2)[1:stride:end]) do i
    idx[] = i
end
