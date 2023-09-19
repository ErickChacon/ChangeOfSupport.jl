# # Continuous inference with observed data at regular grids

# ## Load main packages

using Revise
using ChangeOfSupport
using Meshes
using LinearAlgebra
using SparseArrays
using Statistics
import CairoMakie as MK
import Random

# ## Define a continuous process

Random.seed!(10)

# First, we define the regular B-splines with `q` degree of freedom.
q = 30; order = 3;
bs = RegularBsplines(-100, 100, q, order)

# Evaluate the basis functions
t = range(-100, 100, length = 1000)
B = basis(t, bs)

# The weight of the basis functions are defined by the igmrf.
δ_var = CGMRF(CartesianGrid((q,)), 1, 0.01, 1)
# δ_var = RGMRF(CartesianGrid((q,)), 1, 0.01, 1)
δ = rand(δ_var)

# Now, we can obtain a realization of the continuous stochastic process.
v = B * δ
#-
fig = MK.lines(t, v, label = "W(s)", axis = (title =  "Continuous process",))
MK.scatter!(knotmarks(bs), δ, label = "δ", color = :red)
MK.axislegend()
MK.save("makie.pdf", fig) #src

# ## Generate data observed at aggregated-level

# Lets consider a low resolution regular grid.
tgrid = CartesianGrid(Point(-100.0), Point(100.0), dims = (8,))
Bc = basis(tgrid, bs)
vc = Bc * δ
error = 0.1 * randn(length(vc))
vce = vc + error
#-
fig = MK.rangebars(tgrid, vce, direction = :x, whiskerwidth = 10,
    axis = (title =  "Observed data", xlabel = "s", ylabel = "y"))
MK.save("makie.pdf", fig) #src

# ## Inference

# Inference ignoring the support:

time = centroids(tgrid)[1] |> collect
Bw = basis(time, bs)
P = structure(δ_var)

@time β0, σ²0, κ0 = sample_gam_sparse(vce, Bw, P, 0.1^2, 1, 5000, fix_hyper = false)

# Inference taking into account the support:

@time β, σ², κ = sample_gam_sparse(vce, Bc, P, 0.1^2, 1, 5000, fix_hyper = false)

# Verify MCMC samples.

fig = MK.Figure()
MK.Axis(fig[1,1], title = "β")
MK.series!(β[:, 500:end], solid_color = :gray, linewidth = 0.2)
MK.hlines!(δ, linestyle = :dash, color = :red)
MK.Axis(fig[2,1], title = "κ")
MK.lines!(vec(κ[:, 500:end]), solid_color = :gray)
MK.hlines!([1^2], linestyle = :dash, color = :red)
MK.Axis(fig[3,1], title = "σ²")
MK.lines!(vec(σ²[:, 500:end]), solid_color = :gray)
MK.hlines!([0.1^2], linestyle = :dash, color = :red)
MK.save("makie.pdf", fig) #src

# Comparison with and without ingoring the support.

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
MK.save("makie.pdf", fig) #src
# MK.save("sim-1d-regular.pdf", fig)
