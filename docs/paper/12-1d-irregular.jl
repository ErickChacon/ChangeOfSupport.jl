#' ---
#' title: Change of support in one dimension
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Plots
using Revise
using ChangeOfSupport
using Meshes
using LinearAlgebra
using SparseArrays
using Distributions
using Random

gr(dpi = 250, size = (300*1.5, 200*1.5))

#' ## Define a continuous process
Random.seed!(3)

#' First, we define the regular B-splines with `q` degree of freedom.
q = 50; order = 3;
bs = RegularBsplines(-100, 100, q, order)

#' Evaluate the basis functions
t = range(-100, 100, length = 1000)
B = basis(t, bs)

#' The weight of the basis functions are defined by the igmrf.
δ_var = CGMRF(CartesianGrid((q,)), 1, 0.01, 1)
δ = rand(δ_var)

#' Now, we can obtain a realization of the continuous stochastic process.
v = B * δ

p = plot(t, v, lw = 1, title = " (a) Continuous process", label = "W(s)");
# plot!(igmrf_marks(bs), δ, st = :scatter, color = :red, ms = 2, msw = 0, label = "δ-IGMRF")
savefig(p, "plot.pdf")

#' ## Aggregated-level observation

#' Lets consider a low resolution irregular grid.
Random.seed!(2)
igrid = ChangeOfSupport.RectilinearGrid((-100.0,), (100.0,), (-100 .+ rand(15) * 200,))
Bci = basis(igrid, bs)
vci = Bci * δ
lengths = diff(knotset(igrid)[1])
error = randn(length(vci)) ./ sqrt.(lengths) * 1
vcie = vci + error

function plot_aux(gcenter, gknots, vc, i, p_ms)
    plot(t, v, lw = 1, label = "W(s)", legend = false)
    scatter!(gcenter, vc, xerr = diff(gknots) / 2, ms = p_ms, msc = i, color = i)
    vline!(gknots, c = i, lw = 0.5,ls = :dash)
end

p = plot_aux(centroids(igrid)[1], knotset(igrid)[1], vcie, 4, 4)
savefig(p, "plot.pdf")

#' ## Inference

#' Inference ignoring the support and area

time = centroids(igrid)[1]
Bwi = basis(time, bs)
P = structure(δ_var)

@time β0, σ²0, κ0 = sample_gam(vcie, Matrix(Bwi), Matrix(P), 0.1^2, 1, 3000)

#' Inference ignoring only the support

@time β1, σ²1, κ1 = ChangeOfSupport.sample_gam_area(vcie, Matrix(Bwi), Matrix(P), lengths, 0.1^2, 1, 3000)

#' Inference taking into account the support:

@time β, σ², κ = ChangeOfSupport.sample_gam_area(vcie, Matrix(Bci), Matrix(P), lengths, 0.1^2, 1, 3000)

#' Comparison with and without ingoring the support.

function plot_aux(gcenter, gknots, vc, col, p_ms, β)
    plot(t, B * β[:, 2:end], lw = 0.12, color = :gray, label = false, alpha = 0.2)
    plot!(t, v, lw = 1, label = "W(s)", legend = false, c = 1)
    vline!(gknots, c = :gray, lw = 0.5,ls = :dash, legend = true, label = false)
    scatter!(gcenter, vc[2], xerr = diff(gknots) / 2, ms = p_ms, msc = col[2], color = col[2], label = "observed")
    plot!(t, B * mean(β[:, 2:end], dims = 2), lw = 2, color = col[1],
          label = "predicted mean", alpha = 0.7, legend = true)
end


p1 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
              β[:, 2:3:end]);
title!(p1, "c) Taking into account the support");
p2 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
              β1[:, 2:3:end]);
title!(p2, "b) Ignoring only the support");
p3 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
              β0[:, 2:3:end]);
title!(p3, "a) Ignoring the support and area");
p = plot(p3, p2, p1, layout = (3, 1), size = (500, 900), ylim = (-3, 7))
savefig(p, "plot.pdf")

marks = centroids(igrid)[1]
gridknots = knotset(igrid)[1]
fig = MK.Figure(resolution = (1000, 900))
MK.Axis(fig[1,1], title =  "a) Ignoring the support and homocedasticity")
A = B * β0[:, 100:3:end]
MK.series!(t, A', solid_color = (:gray, 0.1), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1], gridknots[2:end],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β0[:, 100:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(gridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-3, 6)
MK.Axis(fig[1,2], title =  "b) Ignoring the support and heterocedasticity")
A = B * β1[:, 100:3:end]
MK.series!(t, A', solid_color = (:gray, 0.1), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1], gridknots[2:end],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β1[:, 100:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(gridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-3, 6)
MK.Axis(fig[2,1], title =  "c) Considering the support and heterocedasticity")
A = B * β[:, 100:3:end]
[MK.lines!(t, A[:,i], color = (:gray, 0.1), linewidth = 0.3) for i in 1:size(A, 2)];
# MK.series!(t, A', solid_color = (:gray, 0.1), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1], gridknots[2:end],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β[:, 100:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(gridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-3, 6)
MK.axislegend()
MK.ylims!(-2, 8)
# MK.save("sim-1d-irregular.pdf", fig)
MK.save("makie.pdf", fig)


