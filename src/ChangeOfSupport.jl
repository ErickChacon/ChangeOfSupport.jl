module ChangeOfSupport

export RegularKnots, RectilinearGrid, adjacency, difference, structure
using Meshes
import SparseArrays: sparse, spdiagm, I
include("meshes.jl")

# using Random
# using SparseArrays
# using FFTW
# using LinearAlgebra # cholesky
# import LinearAlgebra: (\)
# using Distributions # Gamma
# import SuiteSparse.CHOLMOD: FactorComponent, Dense


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
