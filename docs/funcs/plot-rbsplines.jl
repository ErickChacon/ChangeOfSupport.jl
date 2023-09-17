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

function MK.convert_arguments(P::Type{<:MK.Rangebars}, d::CartesianGrid{1}, x::AbstractVector)
    marks = range(d)[1]
    MK.convert_arguments(P, x, marks[1:end-1], marks[2:end])
end

