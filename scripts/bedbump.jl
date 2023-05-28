include("intro.jl")
using CairoMakie

dx = 200.0
L = 40e3
x = collect(0.0:dx:L)

bedhead = 3e3               # m
bedslope = -0.06            # (m/m)
bedbump_center = 5000.0     # m
bedbump_height = 100.0      # m
bedbump_width = 700         # m

b = bumpy_bed(x, bedhead, bedslope, bedbump_center, bedbump_height, bedbump_width)
lines(b)

