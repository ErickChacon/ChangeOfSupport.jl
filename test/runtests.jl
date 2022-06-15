using Revise
using ChangeOfSupport
using Test
using Meshes
using SparseArrays

const cs = ChangeOfSupport

xxt(x::SparseMatrixCSC) = x .| x'

@testset "RegularKnots" begin
    # regular knots
    knots = RegularKnots(-10, 10, 3, 3, 3)
    @test length(knots) == 11
    @test range(knots) == -25:5:25
    knots = RegularKnots(-10, 10, 3, 0, 0)
    @test range(knots) == -10.0:5.0:10.0

    # get x index
    @test map(x -> cs.get_x_index(x, 1:10), [-1, 1, 9, 11]) == [0, 1, 9, 10]
    @test map(x -> cs.get_x_index(x, 1:10), [2.1, 5.1, 7.9]) == [2, 5, 7]
end

@testset "CartesianGrid" begin
    # 1D - range and centroids
    grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
    @test range(grid) == [-10.0:2.0:10.0]
    @test centroids(grid) == [-9.0:2.0:9.0]
    @test centroidsmat(grid) == reshape(collect(-9.0:2.0:9.0), (10,1))

    # 1D - adjacency
    grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (5,))
    # non cyclic
    @test adjacency(grid, order = 1, cyclic = false) ==
        spdiagm(1 => repeat([true], 4)) |> xxt
    @test adjacency(grid, order = 2, cyclic = false) ==
        spdiagm(2 => repeat([true], 3)) |> xxt
    # cyclic (torus)
    @test adjacency(grid, order = 1, cyclic = true) ==
        spdiagm(1 => repeat([true], 4), 4 => [true]) |> xxt
    @test adjacency(grid, order = 2, cyclic = true) ==
        spdiagm(2 => repeat([true], 3), 3 => repeat([true], 2)) |> xxt

    # 2D - range and centroids
    grid = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (10,20))
    @test range(grid) == [-10.0:2.0:10.0, -10.0:1.0:10.0]
    @test centroids(grid) == [-9.0:2.0:9.0, -9.5:1.0:9.5]
    @test centroidsmat(grid) == [repeat(-9.0:2.0:9.0, 20) repeat(-9.5:1.0:9.5, inner = 10)]

    # 2D - adjacency
    grid = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (3,3))
    # non cyclic
    adjacency(grid, order = 1, cyclic = false) ==
        spdiagm(1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1), 3 => repeat([true], 6)) |> xxt
    adjacency(grid, order = 2, cyclic = false) ==
        spdiagm(2 => sparse([1, 0, 0, 1, 0, 0, 1] .== 1), 6 => repeat([true], 3)) |> xxt
    # cyclic (torus)
    adjacency(grid, order = 1, cyclic = true) ==
        spdiagm(
            1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1),
            2 => sparse([1, 0, 0, 1, 0, 0, 1] .== 1),
            3 => repeat([true], 6),
            6 => repeat([true], 3)) |> xxt
    adjacency(grid, order = 2, cyclic = true) ==
        spdiagm(
            1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1),
            2 => sparse([1, 0, 0, 1, 0, 0, 1] .== 1),
            3 => repeat([true], 6),
            6 => repeat([true], 3)) |> xxt
end

@testset "RectilinearGrid" begin
    # 1D
    ik = [-3.0, 1.0, 8.0, -5.0]
    rgrid = RectilinearGrid((-10.0,), (10.0,), (ik,))
    @test  knotset(rgrid)[1] == [-10.0, -5.0, -3.0, 1.0, 8.0, 10.0]
    @test centroids(rgrid)[1] == [-7.5, -4.0, -1.0, 4.5, 9.0]

    # 2D
    ik1 = [-3.0, 1.0, 8.0, -5.0]
    ik2 = [-7.0, 2.0, 6.0, -2.0]
    rgrid = RectilinearGrid((-10.0,-10.0), (10.0,10.0), (ik1, ik2))
    @test knotset(rgrid)[1] == [-10.0, -5.0, -3.0, 1.0, 8.0, 10.0]
    @test knotset(rgrid)[2] == [-10.0, -7.0, -2.0, 2.0, 6.0, 10.0]
    @test centroids(rgrid)[1] == [-7.5, -4.0, -1.0, 4.5, 9.0]
    @test centroids(rgrid)[2] == [-8.5, -4.5, 0.0, 4.0, 8.0]
end

# @testset "Difference and structure matrices" begin
#     # 1d
#     grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (5,))
#     # non cyclic difference
#     @test difference(grid, order = 1, cyclic = false) ==
#         spdiagm(4, 5, 0 => repeat([-1], 4), 1 => repeat([1], 4))
#     @test difference(grid, order = 2, cyclic = false) ==
#         spdiagm(3, 5, 0 => repeat([1], 3), 1 => repeat([-2], 3), 2 => repeat([1], 3))
#     # cyclic difference
#     @test difference(grid, order = 1, cyclic = true) ==
#         spdiagm(0 => repeat([-1], 5), 1 => repeat([1], 4), -4 => [1])
#     @test difference(grid, order = 2, cyclic = true) ==
#         spdiagm(0 => repeat([1], 5), 1 => repeat([-2], 4), 2 => repeat([1], 3),
#         -3 => repeat([1], 2), -4 => [-2])
#
#     # 2d
#     grid = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (3,3))
#     # non cyclic difference
#     # cyclic difference
# end



