using Revise
using ChangeOfSupport
using Meshes
using Test
using LinearAlgebra

const cs = ChangeOfSupport

@testset "RegularKnots" begin
    knots = RegularKnots(-10, 10, 3, 3, 3)
    @test range(knots) == -25:5:25
    knots = RegularKnots(-10, 10, 3, 0, 0)
    @test range(knots) == -10.0:5.0:10.0
    @test map(x -> cs.get_x_index(x, 1:10), [-1, 1, 9, 11]) == [0, 1, 9, 10]
    @test map(x -> cs.get_x_index(x, 1:10), [2.1, 5.1, 7.9]) == [2, 5, 7]
end

@testset "RegularBsplines" begin
    bs = RegularBsplines(-10, 10, 10, 1)
    @test cs.extendedknots(bs) == RegularKnots(-10.0, 10.0, 9, 0, 1)
    @test cs.boundaryknots(bs) == RegularKnots(-10.0, 10.0, 9, 0, 0)
    @test_throws DomainError(-11, "The values of x should lie in the range of b.") basis(-11, bs)
    @test_throws DomainError(11, "The values of x should lie in the range of b.") basis(11, bs)
    @test basis(-10, bs) == [1.0, zeros(9)...]
    @test basis(-8, bs) == [0.0, 1.0, zeros(8)...]
    @test basis(9.9, bs) == [zeros(9)..., 1.0]
    bs = RegularBsplines(-10, 10, 10, 2)
    @test cs.extendedknots(bs) == RegularKnots(-10.0, 10.0, 8, 1, 2)
    @test cs.boundaryknots(bs) == RegularKnots(-10.0, 10.0, 8, 0, 0)
    @test basis(-10, bs) == [1.0, zeros(9)...]
    @test all(basis(-8, bs)[3:end] .== 0.0)
    @test all(basis(9.9, bs)[1:8] .== 0.0)
    bs = RegularBsplines(-10, 10, 10, 3)
    @test cs.extendedknots(bs) == RegularKnots(-10.0, 10.0, 7, 2, 3)
    @test cs.boundaryknots(bs) == RegularKnots(-10.0, 10.0, 7, 0, 0)
    @test basis(-10, bs) == [0.5, 0.5, zeros(8)...]
    @test all(basis(-8, bs)[4:end] .== 0.0)
    @test all(basis(9.9, bs)[1:7] .== 0.0)
end

@testset "RegularBsplines on vectors" begin
    bs = RegularBsplines(-10, 10, 10, 1)
    @test_throws DomainError([-11, 10], "The values of x should lie in the range of b.") basis([-11, 10], bs)
    @test basis(-10:2:8, bs) == Diagonal(ones(10))
end

@testset "RegularGrid" begin
    grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
    @test range(grid) == [-10.0:2.0:10.0]
    grid = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (10,20))
    @test range(grid) == [-10.0:2.0:10.0, -10.0:1.0:10.0]
end



grid2 = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (10,5))


a = IrregularGrid((1, 1), (10, 10), ([7, 3, 5], [7, 1, 8]))
minimum(a)

a = IrregularGrid((1,), (10,), ([1, 3,],))

a = IrregularGrid2((1,1), (10,10), ([1, 3, 5], [2, 7, 8]))


function IrregularGrid3(minimum::NTuple{Dim,T}, maximum::NTuple{Dim,T},
                        iknots::NTuple{Dim, Vector{T}}) where {Dim,T}
    iknots = map(unique ∘ sort, iknots)
    IrregularGrid3(minimum, maximum, iknots)
end

# a = IrregularGrid3((1,), (10,), ([1, 3, 5],))

struct IrregularGrid4{Dim,T}
    start::Point{Dim,T}
    finish::Point{Dim,T}
    iknots::NTuple{Dim, Vector{T}}
end

# CartesianGrid((100,100),(10.,20.),(1.,1.))

([1, 3, 5, 6, 5], [2, 6, 7, 8, 7]) |> typeof
[[1, 3, 5, 9, 5], [2, 6, 7, 9, 9]] |> typeof

[2, 6, 7, 9, 9] |>
# x -> Point1

Point1[(1,), (2,)]


NTuple{5,NTuple{10}}

(rand(10)..., rand(10)...)



a.minimum
a.maximum
collect(a.iknots...)

CartesianGrid(Point(-10.0), Point(10.0), dims = (40,))

CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T};
              dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, start, (finish - start) ./ dims)

# @testset "RegularBsplines" begin
#     b = RegularBsplines(-10, 10, 3, 6)
#     @test range(allknots(b)) == -20:5:25
#     range(internalknots(b)) == -10:5:10
# end

# Meshes.ivec

# g = CartesianGrid(Point(-10.0), Point(10.0), dims = (13,))
# [coords[1] for coords in coordinates.(vertices(g))] |> diff
#
# g.origin
# d
# CartesianIndices(g.dims .+ 1)
#
#
# ivec(g.origin + (ind.I .- 1) .* g.spacing for ind in inds)
#
#
# g = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (7,13))
# collect(vertices(g))
# g.origin
# g.spacing
# g.dims
# inds = CartesianIndices(g.dims .+ 1)
# (inds[2].I .- 1) .* g.spacing
# Meshes.ivec(1:10)
# Meshes.ivec([1, 4, 5])
#
# function vertices2(g::CartesianGrid)
#   inds = CartesianIndices(g.dims .+ 1)
#   # ivec(g.origin + (ind.I .- 1) .* g.spacing for ind in inds)
#   # [g.origin + (ind.I .- 1) .* g.spacing for ind in inds]
#   knotrange = StepRangeLen(
#     Base.TwicePrecision{Float64}(0),
#     Base.TwicePrecision{Float64}(g.spacing[1]),
#     g.dims[1] + 1)
# end
# Point(1, 0) + vertices2(g)
#
# # g.origin
# # + 
# # (inds[2].I .- 1) .* g.spacing
#
# # [coords[1] for coords in coordinates.(vertices(g))] |> diff
#
#
#
# # grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (7,))
# # knots = [coords[1] for coords in coordinates.(vertices(grid))]
# # diff(knots)
#
# # gridstep = 20.0 / 7
# # knots2 = collect(-10.0:gridstep:10.0)
# # diff(knots2)
# # [knots[i] == knots2[i] for i in 1:7]
#
# # -10.0:7:10.0
# # # -10, 10
#
