module ChangeOfSupport

import Meshes

include("knots.jl")
include("bsplines.jl")
# include("grid.jl")
# include("igmrf.jl")

export RegularKnots, extendedknots, boundaryknots

export RegularBsplines, basis, integral

end # module
