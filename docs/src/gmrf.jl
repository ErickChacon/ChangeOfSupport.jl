#' ---
#' title: Intrinsic Gaussian Markov Random Fields (IGMRF)
#' author: Erick A. Chacón-Montalván
#' weave_options:
#'   term: true
#' ---

using Revise
using Meshes
using Plots
using ChangeOfSupport
using StatsBase
using LinearAlgebra
using Distributions
# using PDMats

gr(dpi = 250)

#' ## One-dimensional IGMRF

#' Defining a grid.
tgrid = CartesianGrid((-10.0,), (10.0,), dims = (5,))
adjacency(tgrid)
structure_cyclic(tgrid)
structure(tgrid)

#' Defining an IGMRF with one neighbour per side.
cgmrf = CGMRF(tgrid, 1, 0.01, 10)
S = structure(cgmrf)
# S[1:5, 1:5]
# structure_base(tgrid, δ = 0.01)
# plot(diag(inv(Matrix(S)))) # vairance 5
gmrf = GMRF(zeros(5000), S, 10)
# X = Matrix{float(eltype(gmrf))}(undef, length(gmrf), 3)
# rand!(gmrf, X)
x = rand(gmrf, (3, 2))

# @time blas = structure(tgrid) + 0.01I
structure_base(tgrid; order = 1)
structure(tgrid)[1:10, 1:10]
structure(tgrid; order = 1)[1:10, 1:10]

sparse(vcat(1:9, 10), vcat(2:10, 1), true, 10, 10)
# mod(1:11, 10)
# (1:10, 1)



# bb = MvNormalCanon(Matrix(gmrf.κ * S))
bb = MvNormalCanon(Matrix(gmrf.κ * S))

sparse(rand(3, 3)) + 1000I
adjacency(tgrid)[1:10, 1:10] + 10I

@time logpdf(gmrf, x)
@time logpdf(bb, x)

bla = cholesky(S)
@time x' * sparse(bla.L);
@time sparse(bla.L)' * x;
@time x' * S * x

@time sum(log.(diag(bla)))
@time logdet(bla) / 2


out = sparse(bla.L) * randn(100)
plot(out[bla.p])


sparse(bla.L)
bla.p



P
= sparse(1:900, bla.p, ones(900))
P{


bla.UP
sparse(bla.PtL)
[1:10, 1:10]

diag(bla)

dense(bla.UP)

@time r = rand(cgmrf, (3, 2));
@time r2 = rand(gmrf, (3, 2));
plot(r)
plot!(r2)
var(r)
var(r2)

# # randn() + randn(rng, reverse(dims)...) * im
a1 = randn(1000) + randn(1000) * im
a2 = randn(ComplexF64, 1000) * sqrt(2)
scatter(a2)
# var(real(a1))
# var(real(a2))


var(r3)
length(gmrf)

rand(gmrf)
# igmrf = IGMRF(tgrid, 1)

x = rand(gmrf)
Distributions._logpdf(gmrf, x)

logpdf(gmrf, x)

typeof(gmrf)
typeof(x)

reshape()
length(grmf)

# Matrix{Vector} <: AbstractMatrix
# size(x, 2)
# structure_base(tgrid)
#
# length(gmrf)
#
# structure_base(sgrid)
# structure(gmrf)

# X = reshape(1:100, 10, 10)
# y = zeros(100)
# copyto!(y, X)

@time bla = structure_base(gmrf.grid; δ = gmrf.δ, order = gmrf.order)
@time bla = structure_base(gmrf.grid)
# @time sparse(bla)

x = zeros(100, 10, 2)
Distributions.rand!(Normal(), x)


plot(fft(bla))

@time x = rand(gmrf, (2, 3))
Matrix{Vector{Float64}}(undef, 2, 3)
[zeros(10) for i in 1:2 for j in 1:3]

@profile x = rand(gmrf, (2, 3))

X = MvNormal(1:3)
x = rand(X, (2, 3))
@eval logpdf(X, x)

size(X) == size(x)
typeof(X)
typeof(x)
Matrix{Vector} <: Matrix


typeof(x)

function plop()
    .x1 = rand(10)
    x1
end

function myfunc()
    A = rand(200, 200, 400)
    maximum(A)
end

@profile myfunc()
Profile.print()

@profile fds = plop()
Profile.print()

plot(x1)
eltype(gmrf)
# x2 = myrand(gmrf)
# plot(x1[:, 2])
# plot!(x2)

bla = randn((2, 3))
[bla...]
vec(bla)

real(3 + 2im) / 3
real((3 + 2im) / 3)
(3 + 2im) / 3


6 / 3 / 2


plot(diag(inv(Matrix(structure(gmrf)))))
# rgmrf = RGMRF(tgrid, 1, 0.00001, 1)
(2 + 0.01) / 10
diag(structure(gmrf))

x = rand(gmrf)
# plot(x)
var(x)
# 1 / (1 + 0.01)

#' Simulate.
x = simulate(igmrf; δ = 0.01, border = 0)

plot(centroids(tgrid)[1], x, lw = 1, title = "IGMRF - N1")

#' Defining an IGMRF with two neighbours per side.
igmrf = IGMRF(tgrid, 2)

#' Simulate.
x = simulate(igmrf; δ = 0.01, border = 0)

plot(centroids(tgrid)[1], x, lw = 1, title = "IGMRF - N2")

#' ## Two-dimensional IGMRF

#' Defining a grid.
sgrid = CartesianGrid((-10.0, -10.0), (10.0, 10.0), dims = (200, 100))

#' Defining an IGMRF with one neighbour per side.
cgmrf = CGMRF(sgrid, 1, 0.01, 1)
x = rand(cgmrf)
x = reshape(x, (100, 200))
heatmap(x)

#' Simulate
x = simulate(igmrf; δ = 0.01, border = 0)

# heatmap(reverse(x, dims = 1), title = "IGMRF - N1")
plot(sgrid, vcat(transpose(reverse(x, dims = 1))...), title = "IGMRF - N1", grid = false,
     size = (500, 518))

#' Defining an IGMRF with 2 neighbours per side.
igmrf = IGMRF(sgrid, 2)

#' Simulate
x = simulate(igmrf; δ = 0.01, border = 0)

plot(sgrid, vcat(transpose(reverse(x, dims = 1))...), title = "IGMRF - N2", grid = false,
     size = (500, 518))

