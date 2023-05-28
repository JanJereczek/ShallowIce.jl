include("intro.jl")
using CairoMakie

N = 51
L = 1.5e6
tspan = (0.0, 20e3)

omega = ComputationDomain(L, N, tspan)
p = Params()
iss = IcesheetState(N)
sstruct = SuperStruct(p, omega, iss)
dt = 10.0

nt = Int(tspan[2] ÷ dt)
ht = zeros(N, nt)
for k in 1:nt
    forwardstep_sia!(sstruct, dt)
    ht[:, k] .= copy(sstruct.iss.h)
end

idx = Observable(1)
stride = 10
hplot = @lift(ht[:, $idx])
fig, ax, p = lines(hplot)
ylims!(ax, (-100, 4000))
record(fig, "sia_N=$N.gif", axes(ht, 2)[1:stride:end]) do i
    idx[] = i
end
# lines(sstruct.iss.h)