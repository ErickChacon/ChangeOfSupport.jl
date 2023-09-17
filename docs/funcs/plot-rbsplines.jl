function boundaryknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, 0, 0)
end

function visualize(bs::RegularBsplines; npoints = 500)
    # evalute bsplines
    t = range(bs.lower, bs.upper, length = npoints)
    Bs = basis(t, bs)
    # plot
    fig = MK.Figure(resolution = (800, 500))
    MK.Axis(fig[1,1], xgridvisible = false)
    MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 1.5);
    [MK.lines!(t, Bs[:,i], label = "$i") for i = 1:size(Bs, 2)]
    MK.current_figure()
end

function visualize(bs::RegularBsplines, t, Bs::AbstractMatrix)
    fig = MK.Figure(resolution = (800, 500))
    MK.Axis(fig[1,1], xgridvisible = false)
    MK.vlines!(range(boundaryknots(bs)), color = :gray, linestyle = :dash, linewidth = 1.5);
    [MK.lines!(t, Bs[:,i], label = "$i") for i = 1:size(Bs, 2)]
    MK.current_figure()
end

function visualize(bs::RegularBsplines, tgrid::CartesianGrid{1}, Bc::AbstractMatrix)
    # evalute bsplines
    t = range(bs.lower, bs.upper, length = 500)
    Bs = basis(t, bs)
    # plot aggregated support
    marks = centroids(tgrid)[1]
    fig = MK.Figure(resolution = (800, 500))
    MK.Axis(fig[1,1], xgridvisible = false)
    for i = 1:size(Bs, 2)
        MK.lines!(t, Bs[:, i])
        MK.rangebars!(Bc[:,i], marks .- step(marks) / 2, marks .+ step(marks) / 2,
            direction = :x, whiskerwidth = 10, linewidth = 2)
    end
    MK.vlines!(range(tgrid)[1], color = :gray, linestyle = :dash, linewidth = 0.8);
    MK.current_figure()
end

function visualize(bs::RegularBsplines, igrid::ChangeOfSupport.RectilinearGrid{1}, Bc::AbstractMatrix)
    # evalute bsplines
    t = range(bs.lower, bs.upper, length = 500)
    Bs = basis(t, bs)
    # plot aggregated support
    gridknots = knotset(igrid)[1]
    # marks = centroids(igrid)[1]
    fig = MK.Figure(resolution = (800, 500))
    MK.Axis(fig[1,1], xgridvisible = false)
    for i = 1:size(Bs, 2)
        MK.lines!(t, Bs[:, i])
        MK.rangebars!(Bc[:,i], gridknots[1:end-1], gridknots[2:end],
            direction = :x, whiskerwidth = 10, linewidth = 2)
    end
    MK.vlines!(gridknots, color = :gray, linestyle = :dash, linewidth = 1.2);
    MK.current_figure()
end

# visualize(bs)
# MK.save("plop.pdf", MK.current_figure())

######################################### RECIPES

# Extension of rangebars to work with CartesianGrid

function MK.convert_arguments(P::Type{<:MK.Rangebars}, d::CartesianGrid{1}, x::AbstractVector)
    marks = range(d)[1]
    MK.convert_arguments(P, x, marks[1:end-1], marks[2:end])
end

# Visualize recipe.

MK.@recipe(Visualize, xpred, ypred, xobs, yobs) do scene
    MK.Attributes(
        labels = nothing
    )
end

MK.MakieLayout.get_plots(plot::Visualize) = plot.plots

function MK.plot!(plot::Visualize)
    xpred = plot[:xpred][]
    ypred = plot[:ypred][]
    xobs = plot[:xobs][]
    yobs = plot[:yobs][]
    [MK.lines!(plot, xpred, ypred[:,i], color = (:gray, 0.1), linewidth = 0.3) for i in 1:size(ypred, 2)];
    # observed
    MK.rangebars!(plot, xobs, yobs, direction = :x, whiskerwidth = 10, linewidth = 2, color = :black);
    MK.scatter!(plot, centroids(xobs)[1], yobs, label = "Observed data", color = :black)
    MK.vlines!(plot, range(xobs)[1], color = :gray, linewidth = 0.5, linestyle = :dash)
    # predictive mean
    MK.lines!(plot, xpred, vec(mean(ypred, dims = 2)) , linewidth = 2, color = :red,
        label = "Predicted mean")
    plot
end


fig = MK.Figure(resolution = (1000, 500))
ax1 = MK.Axis(fig[1,1], title =  "(a) Ignoring the support", xgridvisible = false)
visualize!(t, A0, tgrid, vce)
MK.lines!(t, v, lw = 1)
ax2 = MK.Axis(fig[1,2], title =  "(b) Considering the support", xgridvisible = false)
visualize!(t, A, tgrid, vce)
MK.lines!(t, v, lw = 1, label =  "Latent process")
MK.axislegend()
MK.linkaxes!(ax1, ax2)
MK.hideydecorations!(ax2, grid = false)
MK.ylims!(-2, 8)
MK.save("makie.pdf", fig)

# fig = MK.series(A0[:, 1:10]', solid_color = :gray, linewidth = 0.2, labels = nothing, label = false)
# MK.lines!(rand(1000), linecolor = :red, label = "mean")
# MK.axislegend()
# MK.save("makie.pdf", fig)

# xpred
# ypred
# xobs
# yobs
# # real_process

function bla!(ax::MK.Axis, x, y)
    MK.lines!(ax, x, y, label = "Really?", color = :red, linestyle = :dash)
    # return nothing
end

fig = MK.Figure()
ax = MK.Axis(fig[1,1])
bla!(ax, 1:10, rand(10))
MK.scatter!(1:10, rand(10), label = "common")
MK.axislegend()
MK.save("makie.pdf", fig)
