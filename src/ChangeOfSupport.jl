module ChangeOfSupport

using Meshes

include("knots.jl")
include("meshes.jl")
include("bsplines.jl")
# include("igmrf.jl")

export RegularKnots, extendedknots, boundaryknots

export RegularBsplines, basis, integral, centroids

export IrregularGrid

end # module
