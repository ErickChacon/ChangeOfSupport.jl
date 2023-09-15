# # Regular B-splines

using Revise
using Meshes
using Plots
import CairoMakie as MK
using AlgebraOfGraphics
using ChangeOfSupport

include("funcs/plot-rbsplines.jl")

# gr(dpi = 250)


# ### Define a regular B-splines
#
# Regular B-splines can be defined by providing the boundary knots (minimum and maximum),
# the degree of freedom and de order.

n = 7
bs = RegularBsplines(-10, 10, n, 3)
#-
visualize(bs, npoints = 500)
MK.save("makie.pdf", MK.current_figure()) #src

# ### Evaluating the Bsplines (i)
#
# The basis functions for the Bsplines evaluated at `x` can be obtained with the function
# `basis`.

basis(10.0, bs)

# ### Evaluating the Bsplines (ii)
# Usually we want to evaluate the basis functions at a set of values.

t = range(-10, 10, length = 500)
Bs = basis(t, bs)
#-
visualize(bs, t, Bs)
MK.save("makie.pdf", MK.current_figure()) #src

# ### Evaluating the Bsplines integral from -∞ to x

Bi = integral(t, bs)
#-
visualize(bs, t, Bi)
MK.save("makie.pdf", MK.current_figure()) #src

# When showing the basis functions and the integrals with can see that the integrals reach
# their maximum after the last knot where the basis is not zero.

fig = MK.Figure(resolution = (800, 900))
MK.Axis(fig[1,1], title =  "(a) B-splines")
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 1.5);
[MK.lines!(t, Bs[:,i], label = "$i") for i = 1:size(Bs, 2)]
ax2 = MK.Axis(fig[2,1], title =  "(b) B-splines Integral")
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 1.5);
[MK.lines!(t, Bi[:,i], label = "$i") for i = 1:size(Bi, 2)]
MK.save("makie.pdf", fig) #src

# ### Measure the Bsplines over a regular grid (low resolution)

tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
Bc = basis(tgrid, bs)
# -
visualize(bs, tgrid, Bc)
MK.save("makie.pdf", MK.current_figure()) #src

# ### Measure the Bsplines over a regular grid (high resolution)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (40,))
Bc = basis(tgrid, bs)
# -
visualize(bs, tgrid, Bc)
MK.save("makie.pdf", MK.current_figure()) #src

# ### Measure the Bsplines over a irregular grid (low resolution)
igrid = ChangeOfSupport.RectilinearGrid((-10.0,), (10.0,), (-10 .+ rand(10) * 20,))
Bc = basis(igrid, bs)
# -
visualize(bs, igrid, Bc)
MK.save("makie.pdf", MK.current_figure()) #src

# ### Measure the Bsplines over a irregular grid (high resolution)
igrid = ChangeOfSupport.RectilinearGrid((-10.0,), (10.0,), (-10 .+ rand(40) * 20,))
Bc = basis(igrid, bs)
# -
visualize(bs, igrid, Bc)
MK.save("makie.pdf", MK.current_figure()) #src

