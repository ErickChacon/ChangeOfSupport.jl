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

marks = centroids(bsgrid)[1]
fig = MK.Figure(resolution = (1000, 500))
MK.Axis(fig[1,1], title = "a) Discrete process")
MK.rangebars!(δ, marks .- step(marks) / 2, marks .+ step(marks) / 2,
            direction = :x, whiskerwidth = 10, linewidth = 2)
MK.Axis(fig[1,2], title = "b) Continuous process")
MK.lines!(centroids(tgrid)[1], Bs * δ)
MK.save("makie.pdf", fig)

#' ## Two-dimenional process

n = 50
bsgrid = CartesianGrid((-100.0, -100.0), (100.0, 100.0), dims = (n,n))
bs = NRegularBsplines(bsgrid, 3)
gmrf = RGMRF(bsgrid, 2, 0.1, 1)

tgrid = CartesianGrid((-100.0, -100.0), (100.0, 100.0), dims = (300,200))
Bs = basis(tgrid, bs)
δ = rand(gmrf)

#-

fig = MK.Figure(resolution = (1000, 500))
MK.Axis(fig[1,1], title = "a) Discrete process")
viz!(bsgrid, color = δ, showfacets = true)
MK.Axis(fig[1,2], title = "b) Continuous process")
viz!(tgrid, color = Bs * δ)
MK.save("makie.pdf", fig)

