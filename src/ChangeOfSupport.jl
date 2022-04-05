module ChangeOfSupport

using Meshes
import Meshes
import SparseArrays: sparse, spdiagm
include("meshes.jl")
export RegularKnots, RectilinearGrid, adjacency, difference, structure
export centroids, knotset, centroidsmat

import SuiteSparse # this is only for fixes
import FFTW
import Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra
include("gmrf.jl")
export GMRF, RGMRF, CGMRF

include("bsplines.jl")
export RegularBsplines, NRegularBsplines, basis, integral, boundaryknots, extendedknots
export startingknots, centroidknots

# using Random
# using SparseArrays
# using FFTW
# using LinearAlgebra # cholesky
# import LinearAlgebra: (\)
# using Distributions # Gamma


# include("knots.jl")
# include("rbsplines.jl")
# include("nrbsplines.jl")
# include("gmrf.jl")
# include("igmrf.jl")
# include("cgmrf.jl")
# include("rgmrf.jl")
# include("mcmc.jl")

# export RegularKnots, extendedknots, boundaryknots, startingknots
#
# export RegularBsplines, NRegularBsplines, basis, integral, centroids, centroidsmat, marks
#
# export IrregularGrid, knotset
#
# export IGMRF, ij_to_k, igmrf_precision_base, simulate, igmrf_precision_1, sample_gam, sample_gam_area
#
# export igmrf_knots, igmrf_marks, igmrf_simulate
#
# export sample_gam, igmrf_precision_1_t
#
# export adjacency, adjacency_cyclic, structure_base, CGMRF, precision, structure, structure_cyclic
#
# export RGMRF
#
# export GMRF, plop, old_rand!

end # module
