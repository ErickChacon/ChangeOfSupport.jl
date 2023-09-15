# # N-dimensional Regular B-splines

using Meshes
using Plots
using Revise
using ChangeOfSupport
gr(dpi = 250)

# ## One-dimenional regular b-splines

# Defining the NRegularBsplines
n = 7
bs = NRegularBsplines((-100.0,), (100.0,), (n,), 3)

# Evaluating the Bsplines (i)
basis((10.0,), bs)

# Evaluating the Bsplines (ii)
t = -100.0:100.0 |> collect
Bs = basis(t, bs)
p = plot(t, Bs, legend = false)
savefig(p, "plot.pdf")

# Measure Bsplines over a regular grid (iii)
tgrid = CartesianGrid((-100.0,), (100.0,), dims = (7,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
p = plot(t, Bs, legend = false, grid = nothing, lw = 0.5);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 3);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)
savefig(p, "plot.pdf")

# Measure Bsplines over a rectilinear grid (iv)
igrid = ChangeOfSupport.RectilinearGrid((-100.0,), (100.0,), (-100 .+ rand(10) * 200,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
p = plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)
savefig(p, "plot.pdf")

# # Two-dimensional regular b-splines

# Defining the NRegularBsplines
n1 = 7; n2 = 5
bs = NRegularBsplines((-100.0, -100.0), (100.0, 100.0), (n1, n2), 3)

# Evaluating the Bsplines (i)
basis((10.0, 10.0), bs)

# Evaluating the Bsplines (ii)
sgrid = CartesianGrid(Point(-100.0, -100.0), Point(100.0, 100.0), dims = (100, 80))
coordins = centroidsmat(sgrid)
Bs = basis(coordins, bs)

fig = MK.Figure()
MK.Axis(fig[1,1])
viz!(sgrid, color = Bs[:,1])
MK.Axis(fig[1,2])
viz!(sgrid, color = Bs[:, 9]),
MK.Axis(fig[2,1])
viz!(sgrid, color = Bs[:, 17]),
MK.Axis(fig[2,2])
viz!(sgrid, color = Bs[:, 18]),
MK.save("makie.pdf", fig)

# Measure Bsplines over a regular grid (iii)
sgrid = CartesianGrid(Point(-100.0, -100.0), Point(100.0, 100.0), dims = (10, 10))
Bc = basis(sgrid, bs)

gridknots = range(sgrid)
fig = MK.Figure()
MK.Axis(fig[1,1])
viz!(sgrid, color = Bc[:, 1])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[1,2])
viz!(sgrid, color = Bc[:, 9])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[2,1])
viz!(sgrid, color = Bc[:, 17])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[2,2])
viz!(sgrid, color = Bc[:, 18])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.save("makie.pdf", fig)

# Measure Bsplines over a rectilinear grid (iv)
# fill issue on Meshes.jl abaout sorting axis values
xs = sort(-100 .+ rand(10) * 200)
ys = sort(-100 .+ rand(10) * 200)
xseq = [-100, xs..., 100]
yseq = [-100, ys..., 100]

sgrid = ChangeOfSupport.RectilinearGrid((-100.0, -100.0), (100.0, 100.0), (xs, ys))
Bc = basis(sgrid, bs)

meshesgrid = Meshes.RectilinearGrid(xseq, yseq)
gridknots = knotset(sgrid)

fig = MK.Figure()
MK.Axis(fig[1,1])
viz!(meshesgrid, color = Bc[:,1])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[1,2])
viz!(meshesgrid, color = Bc[:,9])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[2,1])
viz!(meshesgrid, color = Bc[:,17])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.Axis(fig[2,2])
viz!(meshesgrid, color = Bc[:,18])
MK.vlines!(gridknots[1], color = :gray, linestyle = :dash);
MK.hlines!(gridknots[2], color = :gray, linestyle = :dash);
MK.save("makie.pdf", fig)

# gridknots = knotset(sgrid)
# function plot_aux(i)
#     p = heatmap(gridknots[1], gridknots[2], transpose(reshape(Bc[:, i], 11, 11)),
#             grid = false);
#     vline!(p, gridknots[1], c = :gray, lw = 0.5, ls = :dash, label = false);
#     hline!(p, gridknots[2], c = :gray, lw = 0.5, ls = :dash, label = false);
#     p
# end
#
# plot(
#     plot_aux(1),
#     plot_aux(9),
#     plot_aux(17),
#     plot_aux(18),
#     ylim = (-100, 100), xlim = (-100, 100)
# )
