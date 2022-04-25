#' ---
#' title: Intrinsic Gaussian Markov Random Fields (IGMRF)
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Revise
using Meshes
using Plots
using ChangeOfSupport
using Distributions
# using StatsBase
# using LinearAlgebra
# using Distributions
# using SparseArrays
# using PDMats

gr(dpi = 250)

#' ## One-dimensional GMRF

tgrid = CartesianGrid((-10.0,), (10.0,), dims = (100,))
# S = structure(tgrid, δ = 0.01, order = 1)
gmrf = RGMRF(S, 1)
# gmrf = RGMRF(tgrid, 1, 0.01, 1)

@time bla = rand(gmrf, 10);
@time logpdf(gmrf, bla)
plot(bla, legend = false)

#' ## Two-dimensional GMRF

sgrid = CartesianGrid((-10.0,-10.0), (10.0,10.0), dims = (100,120))
S = structure(sgrid, δ = 0.01, order = 2)
gmrf = GMRF(S, 3)
@time bla = rand(gmrf, 3);
@time logpdf(gmrf, bla)
plot(sgrid, bla[:, 1])

reshape(bla[:,1], 100, 120) |>
    transpose |>
    x -> heatmap(range(-10, 10, 100), range(-10, 10, 120), x, aspect_ratio = :equal)
