using DrWatson
@quickactivate "sia"
using OrdinaryDiffEq

# Here you may include files from the source directory
include(srcdir("domain.jl"))
include(srcdir("params.jl"))
include(srcdir("derivatives.jl"))
include(srcdir("physics.jl"))
include(srcdir("utils.jl"))
