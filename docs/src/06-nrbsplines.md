```@meta
EditURL = "https://github.com/ErickChacon/ChangeOfSupport.jl/blob/main/literate/06-nrbsplines.jl"
```

# N-dimensional Regular B-splines

````@example 06-nrbsplines
using Meshes
using Plots
using Revise
using ChangeOfSupport
gr(dpi = 250)
````

## One-dimenional regular b-splines

Defining the NRegularBsplines

````@example 06-nrbsplines
n = 7
bs = NRegularBsplines((-100.0,), (100.0,), (n,), 3)
````

Evaluating the Bsplines (i)

````@example 06-nrbsplines
basis((10.0,), bs)
````

Evaluating the Bsplines (ii)

````@example 06-nrbsplines
t = -100.0:100.0 |> collect
Bs = basis(t, bs)
plot(t, Bs, legend = false)
````

Measure Bsplines over a regular grid (iii)

````@example 06-nrbsplines
tgrid = CartesianGrid((-100.0,), (100.0,), dims = (7,))
Bc = basis(tgrid, bs)

marks = centroids(tgrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, legend = false, grid = nothing, lw = 0.5);
scatter!(marks, Bc, xerr = step(marks) / 2, legend = false,
        msc = cols, color = cols, ms = 3);
vline!(range(tgrid)[1], c = :gray, lw = 0.5, ls = :dash)
````

Measure Bsplines over a rectilinear grid (iv)

````@example 06-nrbsplines
igrid = RectilinearGrid((-100.0,), (100.0,), (-100 .+ rand(10) * 200,))
Bc = basis(igrid, bs)

gridknots = knotset(igrid)[1]
marks = centroids(igrid)[1]
cols = transpose(1:size(Bc)[2])
plot(t, Bs, lw = 1, legend = false, grid = false);
scatter!(marks, Bc, xerr = diff(gridknots) / 2, legend = false,
        msc = cols, color = cols, ms = 4);
vline!(gridknots, c = :gray, lw = 0.5, ls = :dash)
````

# Two-dimensional regular b-splines

Defining the NRegularBsplines

````@example 06-nrbsplines
n1 = 7; n2 = 5
bs = NRegularBsplines((-100.0, -100.0), (100.0, 100.0), (n1, n2), 3)
````

Evaluating the Bsplines (i)

````@example 06-nrbsplines
basis((10.0, 10.0), bs)
````

Evaluating the Bsplines (ii)

````@example 06-nrbsplines
sgrid = CartesianGrid(Point(-100.0, -100.0), Point(100.0, 100.0), dims = (100, 80))
coordins = centroidsmat(sgrid)
Bs = basis(coordins, bs)
````

heatmap(transpose(reshape(Bs[:, 9], 100, 80)))

````@example 06-nrbsplines
plot(
    plot(sgrid, Bs[:, 1]),
    plot(sgrid, Bs[:, 9]),
    plot(sgrid, Bs[:, 17]),
    plot(sgrid, Bs[:, 18]),
    ylim = (-100, 100), xlim = (-100, 100)
)
````

Measure Bsplines over a regular grid (iii)

````@example 06-nrbsplines
sgrid = CartesianGrid(Point(-100.0, -100.0), Point(100.0, 100.0), dims = (10, 10))
Bc = basis(sgrid, bs)

gridknots = range(sgrid)
function plot_aux(i)
    p = plot(sgrid, Bc[:, i]);
    vline!(p, gridknots[1], c = :gray, lw = 0.5, ls = :dash, label = false);
    hline!(p, gridknots[2], c = :gray, lw = 0.5, ls = :dash, label = false);
    p
end

plot(
    plot_aux(1),
    plot_aux(9),
    plot_aux(17),
    plot_aux(18),
    ylim = (-100, 100), xlim = (-100, 100)
)
````

Measure Bsplines over a rectilinear grid (iv)

````@example 06-nrbsplines
sgrid = RectilinearGrid((-100.0, -100.0), (100.0, 100.0),
                      (-100 .+ rand(10) * 200, -100 .+ rand(10) * 200))
Bc = basis(sgrid, bs)

gridknots = knotset(sgrid)
function plot_aux(i)
    p = heatmap(gridknots[1], gridknots[2], transpose(reshape(Bc[:, i], 11, 11)),
            grid = false);
    vline!(p, gridknots[1], c = :gray, lw = 0.5, ls = :dash, label = false);
    hline!(p, gridknots[2], c = :gray, lw = 0.5, ls = :dash, label = false);
    p
end

plot(
    plot_aux(1),
    plot_aux(9),
    plot_aux(17),
    plot_aux(18),
    ylim = (-100, 100), xlim = (-100, 100)
)
````

