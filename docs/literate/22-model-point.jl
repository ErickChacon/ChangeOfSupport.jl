# # Change of support in one dimension

using Revise
using Plots
using ChangeOfSupport
using Meshes
# using LinearAlgebra
using SparseArrays
# using Distributions
using Random
const COS = ChangeOfSupport

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

plot(t, v, lw = 1, title = " (a) Continuous process", label = "W(s)")
plot!(range(centroids(bs)), δ, st = :scatter, color = :red, ms = 2, msw = 0, label = "δ-IGMRF")


#' ## Point-level Estimation

ts = range(-100, 100, length = 50)
Bs = basis(ts, bs)
vs = Bs * δ
es = 0.3 * randn(length(vs))
vse = vs + es
P = structure(δ_var)

plot(ts, vs, st = :scatter, ms = 1.5, msw = 0)
plot!(ts, vse, st = :scatter, ms = 1.5, msw = 0)

X = sparse(Bs)
nonzeros(X) |> length
size(X) |> prod
σ²=0.3^2
κ=1

# @time β, σ², κ = sample_gam(vse, Bs, P, 0.3^2, 1, 1000)
@time out = sample_gam(vse, Bs, P, 0.3^2, 1, 2)
@time out2 = COS.sample_gam_sparse(vse, Bs, P, 0.3^2, 1, 2)

sparse(out2.L) |> Matrix

histogram(out.U \ ones(50), nbins = 1000)
histogram(out2.U \ ones(50), nbins = 1000)

histogram(out2.UP \ collect(1:50), nbins = 20)

histogram(out.L \ collect(1:50), nbins = 20)
histogram(out2.PtL \ collect(1:50), nbins = 20)

sparse(out2.L)[invperm(C.p),:]

out2' |> sparse

out \ reshape(collect(1.0:50.0), 50, 1)
out2 \ collect(1.0:51.0)
out2 \ reshape(collect(1.0:50.0), 50, 1)
out2

typeof(out)
typeof(out2)
UpperTriangular(out2)

# out.U
# Matrix(out2.UP)
# Matrix(out2.L)

# sum(abs.(Matrix(out2) - out))
sum(abs.(out2 - out))

# mcmc
p1 = plot(transpose(σ²), lw = 0.000001, legend = false)
hline!(p1, [0.3^2], linestyle = :dash)

p2 = plot(transpose(κ), lw = 0.000001, legend = false)
hline!(p2, [1^2])

p3 = plot(transpose(β), lw = 0.000001, legend = false)
hline!(δ, linestyle = :dash, color = :gray, lw = 0.1)

plot(p1, p2, p3, layout = (3, 1), size = (500, 600))

# prediction
plot(t, B * β, lw = 0.000001, color = :gray, label = false, alpha = 0.2)
plot!(t, B * δ, label = "process", color = 1)
plot!(ts, vs, st = :scatter, ms = 1.5, msw = 0, label = "no error", c = 2)
plot!(ts, vse, st = :scatter, ms = 1.5, msw = 0, label = "error", c = 3)
