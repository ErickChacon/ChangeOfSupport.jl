using Revise
using Plots
ENV["GKSwstype"] = "100"
using ChangeOfSupport
using Meshes
using LinearAlgebra
using SparseArrays
using Statistics
import CairoMakie as MK
import Random

gr(dpi = 250, size = (300*1.5, 200*1.5))

#' ## Define a continuous process
Random.seed!(10)

#' First, we define the regular B-splines with `q` degree of freedom.
q = 30; order = 3;
bs = RegularBsplines(-100, 100, q, order)

#' Evaluate the basis functions
t = range(-100, 100, length = 1000)
B = basis(t, bs)

#' The weight of the basis functions are defined by the igmrf.
δ_var = CGMRF(CartesianGrid((q,)), 1, 0.01, 1)
# δ_var = RGMRF(CartesianGrid((q,)), 1, 0.01, 1)
δ = rand(δ_var)

# p = plot(δ)
# savefig(p, "plot.pdf")

#' Now, we can obtain a realization of the continuous stochastic process.
v = B * δ

# p = plot(t, v, lw = 1, title = " (a) Continuous process", label = "W(s)");
# # plot!(range(centroids(bs)), δ, st = :scatter, color = :red, ms = 2, msw = 0, label = "δ-IGMRF");
# savefig(p, "plot.pdf")

fig = MK.Figure()
MK.Axis(fig[1,1], title = " (a) Continuous process")
MK.lines!(t, v)
MK.save("makie.pdf", fig)


#' ## Aggregated-level observation

#' Lets consider a low resolution regular grid.
tgrid = CartesianGrid(Point(-100.0), Point(100.0), dims = (8,))
Bc = basis(tgrid, bs)
vc = Bc * δ
error = 0.1 * randn(length(vc))
vce = vc + error
# plot(vce)

#' ## Inference

#' Inference ignoring the support:

time = centroids(tgrid)[1] |> collect
Bw = basis(time, bs)
P = structure(δ_var)

# Qᵦ = Bw'Bw / 0.1^2 + P
# c1 = cholesky(Qᵦ)
# c2 = cholesky(sparse(Qᵦ))

@time β0, σ²0, κ0 = sample_gam_sparse(vce, Bw, P, 0.1^2, 1, 5000, fix_hyper = false)


#' Inference taking into account the support:

@time β, σ², κ = sample_gam_sparse(vce, Bc, P, 0.1^2, 1, 5000, fix_hyper = false)

# p1 = plot(transpose(σ²), lw = 0.000001, legend = false)
# hline!(p1, [0.1^2], linestyle = :dash)
# p2 = plot(transpose(κ), lw = 0.000001, legend = false)
# hline!(p2, [1^2])
# p3 = plot(transpose(β), lw = 0.000001, legend = false)
# hline!(δ, linestyle = :dash, color = :gray, lw = 0.1)
# plot(p1, p2, p3, layout = (3, 1), size = (500, 600))

#' Comparison with and without ingoring the support.

function plot_aux(gcenter, gknots, vc, col, p_ms, β)
    plot(t, B * β[:, 2:end], lw = 0.12, color = :gray, label = false, alpha = 0.2)
    plot!(t, v, lw = 1, label = "Latent process", legend = false, c = 1)
    vline!(gknots, c = :gray, lw = 0.5,ls = :dash, legend = true, label = false)
    scatter!(gcenter, vc[2], xerr = diff(gknots) / 2, ms = p_ms, msc = col[2], color = col[2], label = "Observed data")
    plot!(t, B * mean(β[:, 2:end], dims = 2), lw = 2, color = col[1],
          label = "Predicted mean", alpha = 0.9, legend = :topleft, legendfontsize = 5)
    ylims!((-2, 8))
end

p1 = plot_aux(centroids(tgrid)[1], range(tgrid)[1], (vc, vce), (:red, :black), 3,
              β[:, 500:3:end]);
title!(p1, "(a) Change of support", titlefontsize = 10)
p2 = plot_aux(centroids(tgrid)[1], range(tgrid)[1], (vc, vce), (:red, :black), 3,
              β0[:, 500:3:end]);
title!(p2, "(b) Classical approach", titlefontsize = 10)
p = plot(p1, p2, layout = (1, 2), size = (500*1.5, 200*1.5));
# p
savefig(p, "plot.pdf")
# savefig(p, "figures/application.png")

marks = centroids(tgrid)[1]
fig = MK.Figure(resolution = (1000, 500))
MK.Axis(fig[1,1], title =  "(a) Ignoring the support")
A0 = B * β0[:, 500:3:end]
[MK.lines!(t, A0[:,i], color = (:gray, 0.1), linewidth = 0.3) for i in 1:size(A0, 2)];
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(tgrid)[1], vce, label = "Observed data", color = :black)
MK.rangebars!(vce, marks .- step(marks) / 2, marks .+ step(marks) / 2,
            direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β0[:, 500:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(range(tgrid)[1], color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-2, 8)
# MK.axislegend()
MK.Axis(fig[1,2], title =  "(b) Considering the support")
A = B * β[:, 500:3:end]
# MK.series!(t, A', solid_color = (:gray, 0.1), linewidth = 0.3)
[MK.lines!(t, A[:,i], color = (:gray, 0.1), linewidth = 0.3) for i in 1:size(A, 2)];
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(tgrid)[1], vce, label = "Observed data", color = :black)
MK.rangebars!(vce, marks .- step(marks) / 2, marks .+ step(marks) / 2,
            direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β[:, 500:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(range(tgrid)[1], color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-2, 8)
MK.axislegend()
MK.save("makie.pdf", fig)
# MK.save("sim-1d-regular.pdf", fig)

