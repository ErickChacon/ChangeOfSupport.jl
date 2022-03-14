#' ---
#' title: Regular B-splines
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Meshes
using Plots
using Revise
using ChangeOfSupport
gr(dpi = 250)

#' ### Define a regular B-splines
#' Regular B-splines can be defined by providing the boundary knots (minimum and maximum),
#' the degree of freedom and de order.
n = 7
bs = RegularBsplines(-10, 10, n, 3)

#' ### Evaluating the Bsplines (i)
#' The basis functions for the Bsplines evaluated at `x` can be obtained with the function
#' `basis`.
basis(10.0, bs)

#' ### Evaluating the Bsplines (ii)
#' Usually we want to evaluate the basis functions at a set of values.
t = range(-10, 10, length = 500)
Bs = basis(t, bs)

p1 = plot(t, Bs, lw = 1, legend = false, title = "(a) B-splines");
plot!(p1, range(boundaryknots(bs)), st = :vline, color = :gray, ls = :dash, lw = 0.5)

#' ### Evaluating the Bsplines integral from -∞ to x
Bi = integral(t, bs)

p2 = plot(t, Bi, lw = 1, legend = false, title = "(b) B-splines Integral");
plot!(p2, range(boundaryknots(bs)), st = :vline, color = :gray, ls = :dash, lw = 0.5)

#' When showing the basis functions and the integrals with can see that the integrals reach
#' their maximum after the last knot where the basis is not zero.
plot(p1, p2, layout = (2, 1), size = (600, 800))

#' ### Measure the Bsplines over a regular grid (low resolution)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)

#' ### Measure the Bsplines over a regular grid (high resolution)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (40,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 2);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)

#' ### Measure the Bsplines over a irregular grid (low resolution)
igrid = IrregularGrid((-10.0,), (10.0,), (-10 .+ rand(10) * 20,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)

#' ### Measure the Bsplines over a irregular grid (high resolution)
igrid = IrregularGrid((-10.0,), (10.0,), (-10 .+ rand(40) * 20,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 2);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)
