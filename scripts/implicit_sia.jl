include("intro.jl")
using CairoMakie

N = 51
L = 1.5e6
tspan = (0.0, 20e3)
dt = 10.0

omega = ComputationDomain(L, N, tspan, dt)
p = Params()
iss = IcesheetState(N)
sstruct = SuperStruct(p, omega, iss)
ht = forward_sia(sstruct)

idx = Observable(1)
stride = 10
hplot = @lift(ht[:, $idx])
fig, ax, p = lines(hplot)
ylims!(ax, (-100, 4000))
record(fig, "sia_N=$N.gif", axes(ht, 2)[1:stride:end]) do i
    idx[] = i
end
# lines(sstruct.iss.h)