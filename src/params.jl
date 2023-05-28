struct Params{T<:AbstractFloat}
    A::T
    calving_constant::T
    floating_fraction::T
    g::T
    isostasy_on::Bool
    n::Int
    rho_ice::T
    rho_water::T
    rho_mantle::T
    sliding_constant::T
    tau::T
    fd::T 
    fs::T
end

function Params(;
    A = 1e-16,
    calving_constant = 27.1,
    floating_fraction = 0.15,
    g = 9.81,
    isostasy_on = false,
    n = 3,
    rho_ice = 917.0,
    rho_water = 1000.0,
    rho_mantle = 3300.0,
    sliding_constant = 5.7e-20, # Sliding parameter [Pa^-3 m^2 s^-1]
    tau = 3000.0,
)
    fd = 2*A / (n+2) * (rho_ice * g)^n 
    fs = 10 * fd
    return Params(A, calving_constant, floating_fraction, g, isostasy_on,
        n, rho_ice, rho_water, rho_mantle, sliding_constant, tau, fd, fs)
end

function linear_bed(x::Vector, bedhead::Real, bedslope::Real)
    return bedhead .+ bedslope .* x
end

function bumpy_bed(
    x::Vector,
    bedhead::Real,
    bedslope::Real,
    bedbump_center::Real,
    bedbump_height::Real,
    bedbump_width::Real,
)
    bump = bedbump_height .* exp.( -((x .- bedbump_center) ./ bedbump_width) .^ 2 )
    return linear_bed(x, bedhead, bedslope) + bump
end

"""
#%%
# 1. Model settings
# see 

# Horizontal grid
XGridSpacing     = 200.0        # [m] ori = 100, 200 Hintereisferner and 100 Nigardsbreen
XGridLength      = 40000.0      # [m]  (Determined automaticaly based on XGridSpacing and input length in case of external, Hintereisferner or Nigardsbreen)

# Experiment setup
BedType          = 'Linear'     # 'Linear', 'Bump', 'External', 'Hintereisferner', or 'Nigardsbreen'
WidthType        = 'Constant'   # 'Constant', 'Basin', 'External', 'Hintereisferner', or 'Nigardsbreen'
SideSlopeType    = 'Constant'   # 'Constant', 'External', or 'Hintereisferner'
InitialGlacier   = 'None'       # 'None', 'External', or 'Hintereisferner'
BedDeform        = 'None'       # 'None' or 'Yes'
Calving          = 'None'       # 'None' or 'Floatation' or 'Waterdepth'
"""