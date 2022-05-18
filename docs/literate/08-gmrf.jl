# # Gaussian Markov Random Fields (IGMRF)

using Revise
using Meshes
using Plots
using ChangeOfSupport
using Distributions
using Statistics
const CS = ChangeOfSupport

gr(dpi = 250)

# ## GMRF

# ### One-dimensional

tgrid = CartesianGrid((-10.0,), (10.0,), dims = (500,))
S = structure(tgrid, δ = 0.01, order = 1)
gmrf = GMRF(S, 3)
x = rand(gmrf, 3)
logpdf(gmrf, x)
plot(x)

# Q = 3S
# Σ = inv(Matrix(Q))
# min(diag(Σ)...)
# var(x, dims = 1)
# aux = MvNormalCanon(Matrix(3S))
# logpdf(aux, x)

# ### Two-dimensional GMRF

sgrid = CartesianGrid((-10.0,-10.0), (10.0,10.0), dims = (100,120))
S = structure(sgrid, δ = 0.01, order = 2)
gmrf = GMRF(S, 3)
x = rand(gmrf, 3);
@time logpdf(gmrf, x)
plot(sgrid, x[:, 1])

reshape(x[:,1], 100, 120) |>
    transpose |>
    x -> heatmap(range(-10, 10, 100), range(-10, 10, 120), x, aspect_ratio = :equal)

# ## RGMRF

# ### One-dimensional

rgmrf = RGMRF(tgrid, 2, 0.01, 1)
x = rand(rgmrf, 2)
logpdf(rgmrf, x)
plot(x)

# ### Two-dimensional

rgmrf = RGMRF(sgrid, 2, 0.005, 1)
x = rand(rgmrf, 2)
logpdf(rgmrf, x)
plot(sgrid, x[:, 2])


# ## CGMRF

# ### One-dimensional

cgmrf = CGMRF(tgrid, 2, 0.01, 1)
x = rand(cgmrf, 3)
logpdf(cgmrf, x)
plot(x)


# ### Two-dimensional

sgrid = CartesianGrid((-10.0,-10.0), (10.0,10.0), dims = (100,120))
cgmrf = CGMRF(sgrid, 2, 0.01, 1)
x = rand(cgmrf, 3)
logpdf(cgmrf, x)
plot(sgrid, x[:, 1])

