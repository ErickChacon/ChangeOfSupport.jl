module COSMakieExt

using Meshes
using ChangeOfSupport
using Makie
import Statistics: mean

import ChangeOfSupport: traceplot, traceplot!

# Extension of rangebars for CartesianGrid{1}

function Makie.convert_arguments(P::Type{<:Makie.Rangebars}, d::CartesianGrid{1}, x::AbstractVector)
    marks = range(d)[1]
    Makie.convert_arguments(P, x, marks[1:end-1], marks[2:end])
end

# Extension of rangebars for RectilinearGrid{1}

function Makie.convert_arguments(P::Type{<:Makie.Rangebars}, d::RectilinearGrid{1}, x::AbstractVector)
    marks = knotset(d)[1]
    Makie.convert_arguments(P, x, marks[1:end-1], marks[2:end])
end

# Extension of rangebars for RectilinearGrid{1}

function Makie.convert_arguments(P::Type{<:Makie.Rangebars}, d::GeometrySet{1, Float64, Segment{1, Float64}}, x::AbstractVector)
    marksl = map(s -> coordinates(minimum(s))[1], d)
    marksu = map(s -> coordinates(maximum(s))[1], d)
    Makie.convert_arguments(P, x, marksl, marksu)
end

# Traceplot recipe for MCMC predictions

Makie.@recipe(Traceplot, xpred, ypred, xobs, yobs) do scene
    Makie.Attributes(
        labels = nothing,
        addlines = true
    )
end

Makie.MakieLayout.get_plots(plot::Traceplot) = plot.plots

function Makie.plot!(plot::Traceplot)
    xpred = plot[:xpred][]
    ypred = plot[:ypred][]
    xobs = plot[:xobs][]
    yobs = plot[:yobs][]
    addlines = plot[:addlines][]
    # predictive samples
    [Makie.lines!(plot, xpred, ypred[:,i], color = (:gray, 0.1), linewidth = 0.3) for i in 1:size(ypred, 2)];
    # observed data
    Makie.rangebars!(plot, xobs, yobs, direction = :x, whiskerwidth = 10, linewidth = 2, color = :black, label = "Observed data");
    Makie.scatter!(plot, centroids(xobs)[1], yobs, label = "Observed data", color = :black)
    if addlines
        Makie.vlines!(plot, knotset(xobs)[1], color = :gray, linewidth = 0.5, linestyle = :dash)
    end
    # predictive mean
    Makie.lines!(plot, xpred, vec(mean(ypred, dims = 2)) , linewidth = 2, color = :red, label = "Predicted mean")
    plot
end

end
