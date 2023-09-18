module ChangeOfSupport

using Meshes
import Meshes
import SparseArrays: sparse, sparsevec, spdiagm, spzeros
include("meshes.jl")
export RegularKnots, adjacency, difference, structure
export centroids, knotset, centroidsmat

import SuiteSparse # this is only for fixes
import FFTW
import Distributions: InverseGamma, Gamma, Distributions
import Random: AbstractRNG, randn!
import LinearAlgebra: cholesky, ldiv!, I, LinearAlgebra, Diagonal
include("gmrf.jl")
export GMRF, RGMRF, CGMRF

include("bsplines.jl")
export RegularBsplines, NRegularBsplines, basis, integral, boundaryknots, extendedknots
export startingknots, centroidknots

include("mcmc.jl")
export sample_gam, sample_gam_sparse, sample_gam_area

end
