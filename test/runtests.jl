using ChangeOfSupport
using Test
using Meshes
using SparseArrays
using LinearAlgebra

const cs = ChangeOfSupport

xxt(x::SparseMatrixCSC) = x .| x'

@testset "RegularKnots" begin
    knots = RegularKnots(-10, 10, 3, 3, 3)
    @test length(knots) == 11
    @test range(knots) == -25:5:25
    knots = RegularKnots(-10, 10, 3, 0, 0)
    @test range(knots) == -10.0:5.0:10.0
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
    rgrid = cs.RectilinearGrid((-10.0,), (10.0,), (ik,))
    @test  knotset(rgrid)[1] == [-10.0, -5.0, -3.0, 1.0, 8.0, 10.0]
    @test centroids(rgrid)[1] == [-7.5, -4.0, -1.0, 4.5, 9.0]

    # 2D
    ik1 = [-3.0, 1.0, 8.0, -5.0]
    ik2 = [-7.0, 2.0, 6.0, -2.0]
    rgrid = cs.RectilinearGrid((-10.0,-10.0), (10.0,10.0), (ik1, ik2))
    @test knotset(rgrid)[1] == [-10.0, -5.0, -3.0, 1.0, 8.0, 10.0]
    @test knotset(rgrid)[2] == [-10.0, -7.0, -2.0, 2.0, 6.0, 10.0]
    @test centroids(rgrid)[1] == [-7.5, -4.0, -1.0, 4.5, 9.0]
    @test centroids(rgrid)[2] == [-8.5, -4.5, 0.0, 4.0, 8.0]
end

@testset "Auxiliary variables for RegularBsplines" begin
    # get x index
    @test map(x -> cs.get_index(x, 1:10), [-1, 1, 9, 11]) == [0, 1, 9, 10]
    @test map(x -> cs.get_index(x, 1:10), [2.1, 5.1, 7.9]) == [2, 5, 7]

    # extendedknots
    knots = extendedknots(RegularBsplines(-10.0, 8.0, 10, 2))
    @test range(knots) == -12.0:2.0:12.0
    knots = extendedknots(RegularBsplines(-10.0, 10.0, 10, 3))
    @test range(knots) == -15.0:2.5:17.5
    knots = extendedknots(RegularBsplines(0.0, 14.0, 10, 4))
    @test range(knots) == -6.0:2.0:22.0
end

@testset "RegularBsplines evaluated at scalars" begin
    # ORDER 1
    bs = RegularBsplines(-10.0, 10.0, 10, 1)
    # adequate domain
    errormessage =  "The values of x should lie in the range of bs."
    @test_throws DomainError(-10.1, errormessage) basis(-10.1, bs)
    @test_throws DomainError(10.1, errormessage) basis(10.1, bs)
    # basis at knots
    @test map(x -> basis(x, bs), -10.0:2.0:8.0) ==
        map(x -> SparseVector(10, [x], [1.0]), 1:10)
    @test basis(10.0, bs) == spzeros(10)
    # basis at non-knots
    @test map(x -> basis(x, bs), -9.0:2.0:9.0) ==
        map(x -> SparseVector(10, [x], [1.0]), 1:10)

    # ORDER 2
    bs = RegularBsplines(-10.0, 10.0, 10, 2)
    @test isapprox(basis(-9.0, bs), SparseVector(10, [1, 2], [0.55, 0.45]))
    @test isapprox(basis(-7.0, bs), SparseVector(10, [2, 3], [0.65, 0.35]))
    @test isapprox(basis(-5.0, bs), SparseVector(10, [3, 4], [0.75, 0.25]))
    @test isapprox(basis(5.0, bs), SparseVector(10, [7, 8], [0.25, 0.75]))
    @test isapprox(basis(7.0, bs), SparseVector(10, [8, 9], [0.35, 0.65]))
    @test isapprox(basis(9.0, bs), SparseVector(10, [9, 10], [0.45, 0.55]))

    # ODER 3
    bs = RegularBsplines(-10.0, 10.0, 10, 3)
    @test isapprox(basis(-4.0, bs), SparseVector(10, [3, 4, 5], [0.18, 0.74, 0.08]))
    @test isapprox(basis(-2.0, bs), SparseVector(10, [4, 5, 6], [0.32, 0.66, 0.02]))
    @test isapprox(basis(2.0, bs), SparseVector(10, [5, 6, 7], [0.02, 0.66, 0.32]))
    @test isapprox(basis(4.0, bs), SparseVector(10, [6, 7, 8], [0.08, 0.74, 0.18]))
end

@testset "RegularBsplines evaluated at vectors" begin
    # ORDER 1
    bs = RegularBsplines(-10.0, 10.0, 10, 1)
    # adequate domain
    errormessage =  "The values of x should lie in the range of bs."
    @test_throws DomainError([-10.1, 0.0], errormessage) basis([-10.1, 0.0], bs)
    @test_throws DomainError([10.1, 0.0], errormessage) basis([10.1, 0.0], bs)
    @test_throws DomainError([-10.1, 10.1], errormessage) basis([-10.1, 10.1], bs)
    # basis at knots
    @test basis(-10.0:10.0, bs) ==
        sparse(1:20, repeat(1:10, inner = 2), repeat([1.0], 20), 21, 10)

    # ORDER 2
    bs = RegularBsplines(-10.0, 10.0, 10, 2)
    B = basis(-10.0:2.0:10.0, bs)
    @test isapprox(B, spdiagm(11, 10, 0 => 1.0:-0.1:0.1, -1 => 0.1:0.1:1.0, 1 => [0.0]))

    # ORDER 3
    bs = RegularBsplines(-10.0, 10.0, 10, 3)
    knotstep = step(range(extendedknots(bs)))
    @test basis(-10.0:knotstep:10.0, bs) ≈
        spdiagm(9, 10, 0 => fill(0.5, 9), 1 => fill(0.5, 9), 2 => fill(0.0, 8))
end

@testset "RegularBsplines integrals" begin
    # ORDER 1
    bs = RegularBsplines(-10.0, 10.0, 10, 1)
    @test integral([-10.0], bs) == spzeros(1, bs.df)
    @test integral([10.0], bs) == sparse(fill(1, bs.df), 1:bs.df, fill(step(bs), bs.df))
    @test integral(-9:2:9, bs) == LowerTriangular(fill(step(bs), (10, bs.df)) - LinearAlgebra.I)

    # ORDER 2
    bs = RegularBsplines(-10.0, 10.0, 11, 2)
    @test integral(-10:2:10, bs) == LowerTriangular(fill(step(bs), (11, bs.df)) - LinearAlgebra.I)
end

@testset "RegularBsplines CartesianGrid" begin
    # ORDER 1
    bs = RegularBsplines(-10.0, 10.0, 10, 1)
    tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
    @test basis(tgrid, bs) == Diagonal(fill(1.0, 10))

    # ORDER 2
    bs = RegularBsplines(-10.0, 10.0, 11, 2)
    tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (10,))
    @test basis(tgrid, bs) == spdiagm(0 => fill(0.5, 10), 1 => fill(0.5, 10))[1:10,:]
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

