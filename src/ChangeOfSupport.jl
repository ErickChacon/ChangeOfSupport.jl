module ChangeOfSupport

using Meshes
using SparseArrays
using FFTW

include("knots.jl")
include("meshes.jl")
include("rbsplines.jl")
include("nrbsplines.jl")
include("igmrf.jl")

export RegularKnots, extendedknots, boundaryknots

export RegularBsplines, NRegularBsplines, basis, integral, centroids, centroidsmat

export IrregularGrid, knotset

export IGMRF, ij_to_k, igmrf_precision_base, simulate

end # module
