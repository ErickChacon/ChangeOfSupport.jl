using Revise
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

#' Now, we can obtain a realization of the continuous stochastic process.
v = B * δ

# using Plots
# ENV["GKSwstype"] = "100"
# p = plot(t, v, lw = 1, title = " (a) Continuous process", label = "W(s)");
# # plot!(range(centroids(bs)), δ, st = :scatter, color = :red, ms = 2, msw = 0, label = "δ-IGMRF");
# savefig(p, "plot.pdf")

fig = MK.Figure()
MK.Axis(fig[1,1], title = " (a) Continuous process")
MK.lines!(t, v)
MK.save("makie.pdf", fig) #src

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

A0 = B * β0[:, 500:3:end]
A = B * β[:, 500:3:end]
fig = MK.Figure(resolution = (1000, 500))
ax1 = MK.Axis(fig[1,1], title =  "(a) Ignoring the support", xgridvisible = false)
traceplot!(t, A0, tgrid, vce)
MK.lines!(t, v, lw = 1)
ax2 = MK.Axis(fig[1,2], title =  "(b) Considering the support", xgridvisible = false)
traceplot!(t, A, tgrid, vce)
MK.lines!(t, v, lw = 1, label =  "Latent process")
MK.axislegend()
MK.linkaxes!(ax1, ax2)
MK.hideydecorations!(ax2, grid = false)
MK.ylims!(-2, 8)
MK.save("makie.pdf", fig)
# MK.save("sim-1d-regular.pdf", fig)
