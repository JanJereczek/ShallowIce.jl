
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

function constant_accumulation(t::Real, x::Vector; a::Real = 0.3)
    return fill(a, length(x))
end

function linear_accumulation(t::Real, x::Vector; a1 = 0.25, a2 = -0.4)
    m = (a2-a1) / maximum(x)
    return m .* x .+ a1
end
