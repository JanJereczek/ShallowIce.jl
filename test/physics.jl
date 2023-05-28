include(srcdir("domain.jl"))
include(srcdir("params.jl"))
include(srcdir("derivatives.jl"))
include(srcdir("physics.jl"))

L = 1e6
N = 100
tspan = (0.0, 1e3)

omega = ComputationDomain(L, N, tspan)
p = Params()
iss = IcesheetState(N)
sstruct = SuperStruct(p, omega, iss)

dhdt = fill(0.0, N)
h = vcat( 0.0, fill(1e3, N-2), 0.0 )
forwardstep_sia!(dhdt, h, sstruct, 0.0)