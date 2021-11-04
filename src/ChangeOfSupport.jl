module ChangeOfSupport

using Meshes

include("knots.jl")
include("meshes.jl")
include("bsplines.jl")
include("bsplines-surface.jl")
# include("igmrf.jl")

export RegularKnots, extendedknots, boundaryknots

export RegularBsplines, basis, integral, centroids

export centroidsmat

export NRegularBsplines

end # module
