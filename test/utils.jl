x = collect(0:1e3:1e5)
h = 1e3 .* ( (x .> 3e4) .& (x .< 7e4) )
fig, ax, p = lines(h)
lines!(ax, gaussian_filter(x, H))