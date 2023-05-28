struct Params{T<:AbstractFloat}
    A::T
    rho_ice::T
    g::T
    n::Int
    fd::T 
    fs::T
end

function Params(;
    A = 1e-16,
    rho_ice = 900.0,
    g = 9.81,
    n = 3,
)

    fd = 2*A / (n+2) * (rho_ice * g)^n 
    fs = 10 * fd
    return Params(A, rho_ice, g, n, fd, fs)
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

SecondsPerYear   = 31536000

g                = 9.81         # gracitational acceleration [m s^-2]
        
DensityOfIce     = 917.0        # density of ice [kg m^-3]
DensityOfWater   = 1000.0       # density of water [kg m^-3]
DensityOfRock    = 3300.0       # density of rock [kg m^-3]
        
SlidingParam     = 5.7e-20      # Sliding parameter [Pa^-3 m^2 s^-1]
 
# glacial isostatic rebound       
DeformationParam = 1.9e-24      # Deformation parameter [Pa^-3 s^-1]
TimeDeformation  = 3000.0       # time scale of deformation [years]

# Calving
FloatFrac        = 0.15         # fraction above waterlevel
Proportionality  = 27.1         # Constant of proportionality Calving (Brown, 1982) 


#%%
# 1. Model settings
# see 

# Horizontal grid
XGridSpacing     = 200.0        # [m] ori = 100, 200 Hintereisferner and 100 Nigardsbreen
XGridLength      = 40000.0      # [m]  (Determined automaticaly based on XGridSpacing and input length in case of external, Hintereisferner or Nigardsbreen)

# Time axis
TimeStep         = 0.02         # [years] (its reciprocal must be a whole number)
TotalTime        = 1500.0       # [years]
 
# Experiment setup
BedType          = 'Linear'     # 'Linear', 'Bump', 'External', 'Hintereisferner', or 'Nigardsbreen'
WidthType        = 'Constant'   # 'Constant', 'Basin', 'External', 'Hintereisferner', or 'Nigardsbreen'
SideSlopeType    = 'Constant'   # 'Constant', 'External', or 'Hintereisferner'
InitialGlacier   = 'None'       # 'None', 'External', or 'Hintereisferner'
BedDeform        = 'None'       # 'None' or 'Yes'
Calving          = 'None'       # 'None' or 'Floatation' or 'Waterdepth'

"""