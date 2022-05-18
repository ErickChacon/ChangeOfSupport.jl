# # Meshes functionality in ChangeOfSupport package

using Revise
using Meshes
using ChangeOfSupport
const CS = ChangeOfSupport
using Plots
using LinearAlgebra

# ## Regular Notes

knots = RegularKnots(0, 10, 5, 1, 1)
range(knots) |> collect

CS.get_x_index(3.0, range(knots))

knots2 = RegularKnots(0, 10, 5, 3, 1)
range(knots) |> collect

# ## RectilinearGrid

rg = RectilinearGrid((-10.0,), (10.0,), (-10 .+ rand(5) * 20,))
CS.knotset(rg)[1]
CS.centroids(rg)[1]

rg2 = RectilinearGrid(
    (-10.0, -10.0),
    (10.0, 10.0),
    (-10 .+ rand(5) * 20, -10 .+ rand(5) * 20)
)
CS.knotset(rg2)[1]
CS.knotset(rg2)[2]
CS.centroids(rg2)[1]
CS.centroids(rg2)[2]

# ## CartesianGrid 1d

tgrid = CartesianGrid(Point(-100.0), Point(100.0), dims = (10,))
range(tgrid)
CS.centroids(tgrid)
CS.centroidsmat(tgrid)
adjacency(tgrid)
adjacency(tgrid, order = 2, cyclic = true)

difference(tgrid)
difference(tgrid, order = 2, cyclic = true)

structure(tgrid)
structure(tgrid, δ = 0.001, order = 2, cyclic = true)

# ## CartesianGrid 2d

sgrid = CartesianGrid(Point(-100.0, -100.0), Point(100.0, 100.0), dims = (5,5))
range(sgrid)
CS.centroids(sgrid)
CS.centroidsmat(sgrid)

adjacency(sgrid, order = 1)
difference(sgrid, order = 2, cyclic = false)
structure(sgrid, order = 2, cyclic = false)

