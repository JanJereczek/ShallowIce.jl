function grad()
end

function grad_staggered()
end

function fdx(x::Vector, dx::AbstractFloat, N::Int)
    return (view(x, 2:N) - view(x, 1:N-1)) ./ dx
end

function right_fdx(x::Vector, dx::AbstractFloat, N::Int)
    return (x[2] - x[1]) / dx
end

function central_fdx(x::Vector, dx::AbstractFloat, N::Int)
    return (view(x, 3:N) - view(x, 1:N-2)) ./ (2*dx)
end

function left_fdx(x::Vector, dx::AbstractFloat, N::Int)
    return (x[N] - x[N-1]) / dx
end

function mixed_fdx(x::Vector, dx::AbstractFloat, N::Int)
    return vcat( right_fdx(x, dx, N),
        central_fdx(x, dx, N), left_fdx(x, dx, N) )
end

function divergence(x::Vector, dx::AbstractFloat, N::Int)
    return mixed_fdx(x, dx, N)
end