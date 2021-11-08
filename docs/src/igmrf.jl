#' ---
#' title: Intrinsic Gaussian Markov Random Fields (IGMRF)
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Meshes
using Plots
using Revise
using ChangeOfSupport

gr(dpi = 250)

#' ## One-dimensional IGMRF

#' Defining a grid.
tgrid = CartesianGrid((-10.0,), (10.0,), dims = (200,))

#' Defining an IGMRF with one neighbour per side.
igmrf = IGMRF(tgrid, 1)

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
sgrid = CartesianGrid((-10.0, -10.0), (10.0, 10.0), dims = (200, 200))

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

