#' ---
#' title: Intrinsic Gaussian Markov Random Fields (IGMRF)
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Revise
using Meshes
# using Plots
using ChangeOfSupport
using StatsBase
using LinearAlgebra
using Distributions

gr(dpi = 250)

#' ## One-dimensional IGMRF

#' Defining a grid.
tgrid = CartesianGrid((-10.0,), (10.0,), dims = (1000,))

#' Defining an IGMRF with one neighbour per side.
# igmrf = IGMRF(tgrid, 1)
gmrf = CGMRF(tgrid, 2, 0.01, 1)
structure_base(tgrid)

length(gmrf)

structure_base(sgrid)
structure(gmrf)

x1 = rand(gmrf)
# x2 = myrand(gmrf)
plot(x1)
# plot!(x2)

bla = randn((2, 3))
[bla...]
vec(bla)

real(3 + 2im) / 3
real((3 + 2im) / 3)
(3 + 2im) / 3


6 / 3 / 2


plot(diag(inv(Matrix(structure(gmrf)))))
# rgmrf = RGMRF(tgrid, 1, 0.00001, 1)
(2 + 0.01) / 10
diag(structure(gmrf))

x = rand(gmrf)
# plot(x)
var(x)
# 1 / (1 + 0.01)

#' Simulate.
x = simulate(igmrf; δ = 0.01, border = 0)

plot(centroids(tgrid)[1], x, lw = 1, title = "IGMRF - N1")

#' Defining an IGMRF with two neighbours per side.
igmrf = IGMRF(tgrid, 2)

#' Simulate.
x = simulate(igmrf; δ = 0.01, border = 0)

plot(centroids(tgrid)[1], x, lw = 1, title = "IGMRF - N2")

#' ## Two-dimensional IGMRF

#' Defining a grid.
sgrid = CartesianGrid((-10.0, -10.0), (10.0, 10.0), dims = (200, 100))
randn(3, 4)

#' Defining an IGMRF with one neighbour per side.
igmrf = IGMRF(sgrid, 1)

#' Simulate
x = simulate(igmrf; δ = 0.01, border = 0)

# heatmap(reverse(x, dims = 1), title = "IGMRF - N1")
plot(sgrid, vcat(transpose(reverse(x, dims = 1))...), title = "IGMRF - N1", grid = false,
     size = (500, 518))

#' Defining an IGMRF with 2 neighbours per side.
igmrf = IGMRF(sgrid, 2)

#' Simulate
x = simulate(igmrf; δ = 0.01, border = 0)

plot(sgrid, vcat(transpose(reverse(x, dims = 1))...), title = "IGMRF - N2", grid = false,
     size = (500, 518))

