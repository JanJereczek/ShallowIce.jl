
dx = 0.1
x = 0:dx:5
y = sin.(x)
truth = cos.(x)
fd = mixed_fdx(y, dx, length(x))
error = abs.(fd .- truth)