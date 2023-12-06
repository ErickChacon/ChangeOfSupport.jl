module ChangeOfSupport

using Meshes
import Meshes
import SparseArrays: sparse, sparsevec, spdiagm, spzeros, SparseMatrixCSC
include("meshes.jl")
export RegularKnots, adjacency, difference, structure
export centroids, knotset, centroidsmat

using TimerOutputs

# import SuiteSparse # this is only for fixes
# import FFTW
import Distributions: InverseGamma, Gamma, Distributions, truncated, Normal
# import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra, Diagonal
# include("gmrf.jl")
using GMRFs
export GMRF, CGMRF

include("bsplines.jl")
export RegularBsplines, NRegularBsplines, basis, integral, boundaryknots, extendedknots
export startingknots, centroidknots

include("mcmc.jl")
export sample_gam, sample_gam_sparse, sample_gam_area

# makie extension
function traceplot end
function traceplot! end
export traceplot, traceplot!

# temporal
export knotmarks

end
