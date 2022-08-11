using Revise
using ChangeOfSupport
using Test
using Meshes
using SparseArrays

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

@testset "Auxiliary variables for RegularBsplines" begin
    # get x index
    @test map(x -> cs.get_x_index(x, 1:10), [-1, 1, 9, 11]) == [0, 1, 9, 10]
    @test map(x -> cs.get_x_index(x, 1:10), [2.1, 5.1, 7.9]) == [2, 5, 7]

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
    errormessage =  "The values of x should lie in the range of b."
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
    errormessage =  "The values of x should lie in the range of b."
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

# sparsity
bs = RegularBsplines(-10.0, 10.0, 100000, 3)
vals = -10 .+ 20 * rand(200000)
vals = [-10.0, vals..., 10.0]
@time bla1 = cs.basis_sparse(vals, bs);
@time bla2 = basis(vals, bs);
bla1 == bla2

# integral
bs = RegularBsplines(-10.0, 10.0, 15, 2)
vals = -10 .+ 20 * rand(10)
vals = [-10.0, vals..., 10.0]
@time bla = cs.integral(vals, bs);
@time bla2 = cs.integral_old(vals, bs);

# CartesianGrid
bs = RegularBsplines(-10.0, 10.0, 15, 2)
tgrid = CartesianGrid(Point(-10.0), Point(10.0), dims = (5,))
bla1, ind1 = cs._integral(range(tgrid)[1][1:(end-1)], bs)
bla2, ind2 = cs._integral(range(tgrid)[1][2:end], bs)
cs.sparsebasis(bla1, ind1, bs.df)
cs.sparsebasis(bla2, ind2, bs.df)
bla1t = cs.integral_old(range(tgrid)[1][1:(end-1)], bs)
bla2t = cs.integral_old(range(tgrid)[1][2:end], bs)
bla2t - bla1t

cs.sparsebasis(vcat(bla1[1,:]', bla2), [ind1[1], ind2...], bs.df)
cs.sparsebasis(bla2, ind2, bs.df) - cs.sparsebasis(bla1, ind1, bs.df) +
sparsebasis2(bla2t[1,1], ind1, ind2, 2, bs.df) |>
    x -> x ./ step(range(tgrid)[1])

basis(tgrid, bs)
basis_old(tgrid, bs)

collect(1:10)[1:(end-1)]

# Bc = basis(tgrid, bs)

function sparsebasis2(val::Real, indices1::Vector, indices2::Vector, order::Int, df::Int)
    n = length(indices1)
    reps = indices2 - indices1
    # J = (indices1[i] - order + 1):(indices2[i] - order + 1)
    I = vcat([repeat([i], reps[i]) for i in 1:n]...)
    J = vcat([(indices1[i] - order + 1):(indices2[i] - order) for i in 1:n]...)
    # (n, order) = size(basis)
    # n = length(indices1)
    # I = [i for k in 1:(indices2[i]-indices1[i]+1) for i in 1:n]
    # J = [(indices1[i] - order + k) for k in 1:(indices2[i]-indices1[i]+1) for i in 1:n]
    # sparse(I, J, val, n, df + 1)[:, 1:df]
    # sparsevec(J, val)
    sparse(I, J, val, n, df + 1)[:, 1:df]
end

cs.sparsebasis(bla2, ind2, bs.df) - cs.sparsebasis(bla1, ind1, bs.df)
# +
#     sparsebasis2(1.42857, ind1, ind2, 2, bs.df)



vcat([i:i+3 for i in 1:10]...)

sparsebasis2(10.0, ind1, ind2, 10, 2)

sparsevec([1, 3, 4], 4.0)



# @time bla2, ind2 = cs._integral(vals, bs);

# range(extendedknots(bs)) |> collect
# step(range(extendedknots(bs)))
# vals = -10 .+ 20 * rand(5)
vals = [-2]
# vals = [-10.0, vals..., 10.0]
# cs.get_x_index(jj)

@time bla2, ind2 = cs._integral(vals, bs);
@time bla1 = cs.integral(vals, bs);
bla1
cs.sparsebasis(bla2, ind2, bs.df)

cs.nonzerobasis(vals, bs)[1]
# bla1 == bla2
#
Matrix(bla1)





length(range(extendedknots(bs)))
# df + 

# @time bla = findall(x -> x < 2, 1:10000000)

x = 10.0
bla1 = cs.basis_sparse(x, bs)
bla2 = basis(x, bs)
bla1 == bla2

# @time bla2, jx2 = cs._basis(vals, bs);
# cs._sparsebasis(bla2, jx2, bs.df)



ss = rand(3, 4)
cat(ss, dims = 2)

vec(ss)

# bla1
#
# [(i+2) for i in 1:10]


# for i = 1:length(jx2)
# end

Matrix(bla1)
bla2

@time in1 = integral(-10.0:0.9:8.0, bs) |> Matrix;
@time in2 = cs._integral(-10.0:0.9:8.0, bs);

    # knotstep * 
reverse(cumsum(reverse(in2, dims = 2), dims = 2), dims = 2)


# bs2 = RegularBsplines(-10.0, 10.0, 40 + 1, 3 + 1)
# @time ibasis1 = basis(-10.0:0.01:10.0, bs);
# @time ibasis2 = cs.basis_old(-10.0:0.01:10.0, bs);
@time bla = integral(-10.0:2.0:10.0, bs)

# integral(2.0, bs)

# ibasis1 == ibasis2


# [2:bs.df]

@time bla = reverse(ibasis)

@time cumsum(ibasis[end:-1:1])[end:-1:1]
@time reverse(cumsum(reverse(ibasis)))

# bla[1] == 0.0


cumsum(reverse(ibasis))


Matrix(integral(-9.0, bs) .+ 1.0)

    b = 
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis)))





# bla = basis(1.0, bs)
# Matrix(bla)


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

