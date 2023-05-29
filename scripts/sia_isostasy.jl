include("intro.jl")
using CairoMakie

N = 51
L = 1.5e6
tspan = (0.0, 20e3)
dt = 1.0

omega = ComputationDomain(L, N, tspan, dt)
p = Params(isostasy_on = true)
h0, b0 = fill(0.0, N), fill(0.0, N)
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
record(fig, "plots/sia_isostasy_N=$N.gif", axes(ht, 2)[1:stride:end]) do i
    idx[] = i
end
