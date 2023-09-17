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
Random.seed!(7)

#' First, we define the regular B-splines with `q` degree of freedom.
q = 40; order = 3;
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
# Random.seed!(19)
Random.seed!(1)
igrid = ChangeOfSupport.RectilinearGrid((-100.0,), (100.0,), (-100 .+ rand(25) * 200,))

Bci = basis(igrid, bs)
index = sample(axes(Bci, 1), 15, replace = false, ordered = true)
Bci = Bci[index,:]
vci = Bci * δ
lengths = diff(knotset(igrid)[1])[index]
error = randn(length(vci)) ./ sqrt.(lengths) * 1
vcie = vci + error

function plot_aux(gcenter, gknots, vc, i, p_ms, index)
    plot(t, v, lw = 1, label = "W(s)", legend = false)
    scatter!(gcenter[index], vc, xerr = diff(gknots)[index] / 2, ms = p_ms, msc = i, color = i)
    vline!(gknots, c = i, lw = 0.5,ls = :dash)
end

plot_aux(centroids(igrid)[1], knotset(igrid)[1], vcie, 4, 4, index)

#' ## Inference

#' Inference ignoring the support and area

time = centroids(igrid)[1][index]
Bwi = basis(time, bs)
P = structure(δ_var)

@time β0, σ²0, κ0 = sample_gam(vcie, Matrix(Bwi), Matrix(P), 0.1^2, 1, 4000)

#' Inference ignoring only the support

@time β1, σ²1, κ1 = sample_gam_area(vcie, Matrix(Bwi), Matrix(P), lengths, 0.1^2, 1, 4000)

#' Inference taking into account the support:

@time β, σ², κ = sample_gam_area(vcie, Matrix(Bci), Matrix(P), lengths, 0.1^2, 1, 4000)

#' Comparison with and without ingoring the support.

marks = centroids(igrid)[1]
gridknots = knotset(igrid)[1]
ugridknots = unique([gridknots[1:end-1][index]..., gridknots[2:end][index]...])
fig = MK.Figure(resolution = (1000, 900))
MK.Axis(fig[1,1], xgridvisible = false,
    title =  "a) Ignoring the support and homocedasticity")
A = B * β0[:, 500:3:end]
MK.series!(t, A', solid_color = (:gray, 0.2), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1][index], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1][index], gridknots[2:end][index],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β0[:, 500:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(ugridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
MK.ylims!(-8, 3.5)
MK.Axis(fig[1,2], xgridvisible = false,
    title =  "b) Ignoring the support and heterocedasticity")
A = B * β1[:, 500:3:end]
MK.series!(t, A', solid_color = (:gray, 0.2), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1][index], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1][index], gridknots[2:end][index],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β1[:, 500:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(ugridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
# MK.ylims!(0, 6)
# MK.ylims!(-2, 7)
MK.ylims!(-8, 3.5)
MK.Axis(fig[2,1], xgridvisible = false,
    title =  "c) Considering the support and heterocedasticity")
A = B * β[:, 500:3:end]
MK.series!(t, A', solid_color = (:gray, 0.2), linewidth = 0.3)
MK.lines!(t, v, lw = 1, label = "Latent process")
MK.scatter!(centroids(igrid)[1][index], vcie, label = "Observed data", color = :black)
MK.rangebars!(vcie, gridknots[1:end-1][index], gridknots[2:end][index],
    direction = :x, whiskerwidth = 10, linewidth = 2, color = :black)
MK.lines!(t, vec(B * mean(β[:, 500:3:end], dims = 2)) , linewidth = 2, color = :red,
    label = "Predicted mean")
MK.vlines!(ugridknots, color = :gray, linewidth = 0.5, linestyle = :dash)
# MK.ylims!(0, 6)
# MK.ylims!(-2, 7)
MK.ylims!(-8, 3.5)
# # MK.save("sim-1d-irregular.pdf", fig)
MK.save("makie.pdf", fig)



# function plot_aux(gcenter, gknots, vc, col, p_ms, β, index)
#     plot(t, B * β[:, 2:end], lw = 0.12, color = :gray, label = false, alpha = 0.2)
#     plot!(t, v, lw = 1, label = "W(s)", legend = false, c = 1)
#     vline!(gknots, c = :gray, lw = 0.5,ls = :dash, legend = true, label = false)
#     scatter!(gcenter[index], vc[2], xerr = diff(gknots)[index] / 2, ms = p_ms, msc = col[2], color = col[2], label = "observed")
#     plot!(t, B * mean(β[:, 2:end], dims = 2), lw = 2, color = col[1],
#           label = "predicted mean", alpha = 0.8, legend = true)
# end
#
#
# p1 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
#               β[:, 2:3:end], index);
# title!(p1, "c) Taking into account the support");
# p2 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
#               β1[:, 2:3:end], index);
# title!(p2, "b) Ignoring only the support");
# p3 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
#               β0[:, 2:3:end], index);
# title!(p3, "a) Ignoring the support and area");
# plot(p3, p2, p1, layout = (3, 1), size = (500, 900), ylim = (-6, 13))





# # #' ## Aggregated-level Estimation
# # g = CartesianGrid((15))
# # idx = sample(1:nelements(g), 10, replace = false) |> sort
# # g_s = map(x -> element(g, x), idx)
# # # col = repeat(["white"], 15)
# # # col[idx] .= "red"
# # gcenter = vcat(coordinates.(centroid.(g_s))...)
# # gmin = vcat(coordinates.(minimum.(g_s))...)
# # gmax = vcat(coordinates.(maximum.(g_s))...)
# #
# #
# #
# # coordinates(minimum(g_s[1]))[]
# # typeof(g_s)
# #
# # integral(gmax, bs)
# #
# # plot(gcenter, ones(10), xerr = 0.5, st = :scatter)
# #
# # plot(g, c = col)
# #
# #
# # plot(g)
# #
# # sample(CartesianGrid((5)), 4)
# #
# #
# # |>
# #     elements |>
# #     collect
#
# #' Lets consider a low resolution regular grid.
# Random.seed!(3)
# igrid = IrregularGrid((-100.0,), (100.0,), (-100 .+ rand(15) * 200,))
# Bci = basis(igrid, bs)
# vci = Bci * δ
# lengths = diff(knotset(igrid)[1])
# error_sigmas = sqrt(4) ./ sqrt.(lengths)
# error = randn(length(vci)) .* error_sigmas
# # error = randn(length(vci)) * 10
# vcie = vci + error
# scatter(centroids(igrid)[1], vci, xerr = diff(knotset(igrid)[1]) / 2, ms = 2)
# scatter!(centroids(igrid)[1], vcie, xerr = diff(knotset(igrid)[1]) / 2, ms = 2)
# # plot(vci)
# # scatter!(vcie)
#
#
# time = centroids(igrid)[1]
# Bwi = basis(time, bs)
# P = structure(δ_var)
#
# # @time β0, σ²0, κ0 = sample_gam(vcie, Bwi, P, 0.1^2, 1, 3000)
# # @time β, σ², κ = sample_gam(vcie, Bci, P, 0.1^2, 1, 3000)
# aa = sample_gam_area(vcie, Bci, P, lengths, 4, 1, 2)
# heatmap(aa)
# cholesky(aa)
# # inv(aa)
#
# D1 = Diagonal(sqrt.(lengths))
# Bci2 = D1 * Bci
# S = transpose(Bci2) * Bci2 / 4 + P
# heatmap(S)
# cholesky(S)
#
# # bla = [1 0 0; 0 0 0; 0 0 3]
# # inv(bla)
#
# # S = transpose(Bci) * Diagonal(lengths) * Bci * 4 + P
# heatmap(S)
# heatmap(Bci)
# cholesky(S)
#
# cholesky()
#
#
#
# # p1 = plot(transpose(σ²), lw = 0.000001, legend = false)
# # hline!(p1, [0.1^2], linestyle = :dash)
# # p2 = plot(transpose(κ), lw = 0.000001, legend = false)
# # hline!(p2, [1^2])
# # p3 = plot(transpose(β), lw = 0.000001, legend = false)
# # hline!(δ, linestyle = :dash, color = :gray, lw = 0.1)
# # plot(p1, p2, p3, layout = (3, 1), size = (500, 600))
#
# # ylim = [-5, 5]
# function plot_aux(gcenter, gknots, vc, col, p_ms, β)
#     plot(t, B * β[:, 2:end], lw = 0.12, color = :gray, label = false, alpha = 0.2)
#     plot!(t, v, lw = 1, label = "W(s)", legend = false, c = 1)
#     vline!(gknots, c = :gray, lw = 0.5,ls = :dash, legend = true, label = false)
#     scatter!(gcenter, vc[2], xerr = diff(gknots) / 2, ms = p_ms, msc = col[2], color = col[2], label = "observed")
#     plot!(t, B * mean(β[:, 2:end], dims = 2), lw = 2, color = col[1],
#           label = "predicted mean", alpha = 0.6, legend = true)
# end
#
#
# p1 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
#               β[:, 2:3:end]);
# p2 = plot_aux(centroids(igrid)[1], knotset(igrid)[1], (vci, vcie), (:red, :black), 3,
#               β0[:, 2:3:end]);
# # pp = plot(p1, p2, layout = (1, 2), size = (500 * 1.5, 220 * 1.5))
# plot(p1, p2, layout = (2, 1), size = (500, 600), ylim = (-6, 10))
#
# # plot(p1, p2, layout = (2, 1), size = (500, 600), ylim = (-12, 17)) # 7
# # plot(p1, p2, layout = (2, 1), size = (500, 600), ylim = (-10, 15)) # 8
# # plot(p1, p2, layout = (2, 1), size = (500, 600), ylim = (-7, 12)) # 13
#
# # plot(p1, p2, layout = (2, 1), size = (500, 600), ylim = (-5, 15))
#
# # #' ## Aggregated-level Estimation
# #
# # #' Lets consider a low resolution regular grid.
# # tgrid = CartesianGrid(Point(-100.0), Point(100.0), dims = (13,))
# # Bc = basis(tgrid, bs)
# # vc = Bc * δ
# # error = 0.1 * randn(length(vc))
# # vce = vc + error
# #
# # time = centroids(tgrid)[1] |> collect
# # Bw = basis(time, bs)
# # P = igmrf_precision_1_t(q; δ = 0.01)
# # @time β0, σ²0, κ0 = sample_gam(vce, Bw, P, 0.1^2, 1, 3000)
# #
# # @time β, σ², κ = sample_gam(vce, Bc, P, 0.1^2, 1, 3000)
# #
# # p1 = plot_aux(centroids(tgrid)[1], range(tgrid)[1], (vc, vce), (:red, :black), 3,
# #               β[:, 2:3:end]);
# # p2 = plot_aux(centroids(tgrid)[1], range(tgrid)[1], (vc, vce), (:red, :black), 3,
# #               β0[:, 2:3:end]);
# # pp = plot(p1, p2, layout = (1, 2), size = (500 * 1.5, 220 * 1.5))
# #
# # # savefig(pp, "plop.pdf")
# #
# # # plot(p1, p2, layout = (1, 2), size = (900 / 2, 400 /2), dpi = 100)
# #
# # # gr(dpi = 250, size = (300*1.5*2/2, 200*1.5/2))
# # # plot(p1, p2, layout = (1, 2),)
# # # plot(p1, p2, layout = (1, 2), size = (100, 60))
# #
# # # fig = plot(plot(1:10), plot(1:10), layout = (1, 2), size = (500*2, 220*2), dpi = 200)
# # # savefig("plop.pdf")
# #
# # # savefig(p1, "figures/p1.png")
# # # savefig(p2, "figures/p2.png")
# # #
# #
# #
# # # plot(p1, p2, layout = (2, 1), size = (500, 600))
# #
# # # mean(InverseGamma(0.001, 0.001))
# # # mean(1 ./ rand(Gamma(0.001, 1 / 1000), 1000))
# #
# # # Gamma(0.001, 1 / 1000)
# #
# # # heatmap(transpose(Bc) * Bc)
# # # cholesky(transpose(B) * B)
# # # heatmap(transpose(B) * B)
# #
# # # x = 0:0.001:100
# # # # g = Gamma(0.01, 1 / 0.005)
# # # g = Gamma(10.1 * 1000, 10 / 10.1 / 1000)
# # # f = pdf(g, x)
# # # plot(x, f, lw = 0.1)
# # # mean(g)
# #
# # # mean(Gamma(0.001, 1/ 0.001))
# # # mean(Gamma(0.001, 1/ 0.001))
# #
# # # plot_aux(centroids(tgrid)[1], range(tgrid)[1], (vc, vce), (2, 3), 4)
# # # plot!(t, B * β[:, 3500:end], lw = 0.000001, color = :gray, label = false, alpha = 0.4)
# # # plot!(centroids(tgrid)[1], Bc * β[:, 3500:end], color = :red, label = false, st = :scatter, ms = 0.7)
# #
# # # sum((vce - Bc * β[:, 1]).^2) / length(vce)
# #
# # # plot!(centroids(tgrid)[1], Bc * a, color = :red)
# #
# # # plot!(centroids(tgrid)[1], a)
# #
