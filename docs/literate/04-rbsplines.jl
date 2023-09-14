# # Regular B-splines

using Revise
using Meshes
using Plots
import CairoMakie as MK
using AlgebraOfGraphics
using ChangeOfSupport
gr(dpi = 250)

# ### Define a regular B-splines
#
# Regular B-splines can be defined by providing the boundary knots (minimum and maximum),
# the degree of freedom and de order.
n = 7
bs = RegularBsplines(-10, 10, n, 3)

# ### Evaluating the Bsplines (i)
#
# The basis functions for the Bsplines evaluated at `x` can be obtained with the function
# `basis`.

basis(10.0, bs)

# ### Evaluating the Bsplines (ii)
# Usually we want to evaluate the basis functions at a set of values.
t = range(-10, 10, length = 500)
Bs = basis(t, bs)

# Bss = sparse(Bs)
# length(findall(!iszero, Bss)) / prod(size(Bss))

function boundaryknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, 0, 0)
end

ENV["GKSwstype"] = "100"
# GKS: can't connect to GKS socket application
# qt.qpa.xcb: XKeyboard extension not present on the X server
# There is no XRandR 1.2 and later version available. There will be only fake screen(s) to use.
# connect: Connection refused
# GKS: can't connect to GKS socket application

p1 = plot(t, Bs, lw = 1, legend = false, title = "(a) B-splines");
plot!(p1, range(boundaryknots(bs)), st = :vline, color = :gray, ls = :dash, lw = 0.5);
savefig(p1, "plot.pdf")

ggplot_theme = MK.Theme(
    Axis = (
        backgroundcolor = :gray90,
        # backgroundcolor = :white,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
    )
)

simple_theme = MK.Theme(
    Axis = (
        # backgroundcolor = :gray90,
        # # backgroundcolor = :white,
        # leftspinevisible = false,
        # rightspinevisible = false,
        # bottomspinevisible = false,
        # topspinevisible = false,
        # xgridcolor = :white,
        # ygridcolor = :white,
        xgridvisible = false,
    )
)


# MK.set_theme!(ggplot_theme)

# fig = MK.series(t, Bs', color = :Set1, axis = (title =  "(a) B-splines",));
fig = MK.series(t, Bs', axis = (title =  "(a) B-splines",));
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 0.8);
MK.save("makie.pdf", fig) #src

fig = MK.Figure(resolution = (800, 500))
MK.Axis(fig[1,1], title =  "(a) B-splines", xgridvisible = false)
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 1);
[MK.lines!(t, Bs[:,i]) for i = 1:size(Bs, 2)]
MK.save("makie.pdf", fig) #src

using Tables
df = Tables.table(hcat(t, Bs))
plt = data(df) * mapping(1 => "t", 2:8 .=> :value, color = dims(1) => renamer(2:8)) * visual(MK.Lines) +
    data((x = range(boundaryknots(bs)),)) * mapping(:x) * visual(MK.VLines, linestyle = :dash)
fig = draw(plt, figure = (resolution = (800, 500),))
MK.save("makie.pdf", fig) #src



# ### Evaluating the Bsplines integral from -∞ to x
Bi = integral(t, bs)

# p2 = plot(t, Bi, lw = 1, legend = false, title = "(b) B-splines Integral");
# plot!(p2, range(boundaryknots(bs)), st = :vline, color = :gray, ls = :dash, lw = 0.5)

fig = MK.series(t, Bi', color = :Set1, axis = (title =  "(b) B-splines Integral",));
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 0.8);
MK.save("makie.pdf", fig) #src

# When showing the basis functions and the integrals with can see that the integrals reach
# their maximum after the last knot where the basis is not zero.
# p = plot(p1, p2, layout = (2, 1), size = (600, 800));
# savefig(p, "plot.pdf")

fig = MK.Figure()
ax1 = MK.Axis(fig[1,1], title =  "(a) B-splines")
MK.series!(ax1, t, Bs', color = :Set1);
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 0.8);
ax2 = MK.Axis(fig[2,1], title =  "(b) B-splines Integral")
MK.series!(ax2, t, Bi', color = :Set1);
MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 0.8);
MK.save("makie.pdf", fig) #src

# ### Measure the Bsplines over a regular grid (low resolution)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
p = plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)
savefig(p, "plot.pdf")

# fig = MK.series(t, Bs')
fig = MK.Figure()
MK.Axis(fig[1, 1])
# MK.scatter!(marks, Bc[:,1])
for i = 1:size(Bs, 2)
    MK.lines!(t, Bs[:, i])
    MK.rangebars!(Bc[:,i], marks .- step(marks) / 2, marks .+ step(marks) / 2, direction = :x, label = "Basis $i",
        whiskerwidth = 10)
end
MK.vlines!(range(tgrid)[1], color = :gray, linestyle = :dash, linewidth = 0.8);
MK.axislegend()
MK.save("makie.pdf", fig) #src

# using Tables
# plt = data(Tables.table(Bs)) * mapping([:Column1, :Column2, :Columns])
#
# fig = MK.series(t, Bs')
# # MK.scatter!(marks, Bc[:,1])
# # MK.vlines!(range(tgrid)[1])
# MK.save("makie.pdf", fig) #src
#



# ### Measure the Bsplines over a regular grid (high resolution)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (40,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 2);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)

# ### Measure the Bsplines over a irregular grid (low resolution)
igrid = CS.RectilinearGrid((-10.0,), (10.0,), (-10 .+ rand(10) * 20,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)

# ### Measure the Bsplines over a irregular grid (high resolution)
igrid = CS.RectilinearGrid((-10.0,), (10.0,), (-10 .+ rand(40) * 20,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 2);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)
