module ChangeOfSupport

using Meshes

include("knots.jl")
include("meshes.jl")
include("rbsplines.jl")
include("nrbsplines.jl")
# include("igmrf.jl")

export RegularKnots, extendedknots, boundaryknots

export RegularBsplines, NRegularBsplines, basis, integral, centroids, centroidsmat

export IrregularGrid, knotset

end # module
