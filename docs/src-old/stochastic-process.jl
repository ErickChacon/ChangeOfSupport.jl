#' ---
#' title: N-dimensional Regular B-splines
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Meshes
using Plots
using Revise
using ChangeOfSupport
gr(dpi = 250)

#' ## One-dimensional process

n = 30
bsgrid = CartesianGrid((-100.0,), (100.0,), dims = (n,))
bs = NRegularBsplines(bsgrid, 3)
gmrf = RGMRF(bsgrid, 2, 0.01, 1)

tgrid = CartesianGrid((-100.0,), (100.0,), dims = (100,))
Bs = basis(tgrid, bs)
δ = rand(gmrf)
plot(plot(δ), plot(Bs * δ), ylim = (-10, 10))

#' ## Two-dimenional process

n = 50
bsgrid = CartesianGrid((-100.0, -100.0), (100.0, 100.0), dims = (n,n))
bs = NRegularBsplines(bsgrid, 3)
gmrf = RGMRF(bsgrid, 2, 0.1, 1)
tgrid = CartesianGrid((-100.0, -100.0), (100.0, 100.0), dims = (300,200))
Bs = basis(tgrid, bs)
δ = rand(gmrf)
plot(plot(bsgrid, δ), plot(tgrid, Bs * δ))

