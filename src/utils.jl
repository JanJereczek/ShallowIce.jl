
function gaussian_filter(x::Vector, y::Vector)
    
    if length(x) != length(y)
        error("Vector lengths do not match.")
    end

    yfilt = similar(y)
    sigma = 100e3
    for i in eachindex(x)
        yfilt[i] = sum(y .* normed_gaussian(x, x[i], sigma) )
    end
    return yfilt
end

function normed_gaussian(x::Vector, mu::Real, sigma::Real)
    gaussian = exp.( -0.5 .* ( (x .- mu)./sigma ) .^ 2 )
    return gaussian ./ sum(gaussian)
end

function constant_accumulation(sstruct::SuperStruct; a::Real = 0.3)
    return fill(a, sstruct.omega.N)
end

function linear_accumulation(sstruct::SuperStruct; a1 = 0.25, a2 = -0.4)
    m = (a2-a1) / maximum(sstruct.omega.xH)
    return m .* sstruct.omega.xH .+ a1
end

function zero_bc!(x::Vector{T}) where {T<:Real}
    x[0] = T(0)
    x[end] = T(0)
    return nothing
end