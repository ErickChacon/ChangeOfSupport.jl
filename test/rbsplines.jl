using Revise
using ChangeOfSupport
using Meshes
using Test
using BenchmarkTools

const cs = ChangeOfSupport

n = 10
bs = RegularBsplines(-10.0, 10.0, n, 3)
range(extendedknots(bs))
# |> length

coordinates.(vertices(cs.cartesiangrid(bs)))

coordinates(Point((1.,)))


# CartesianGrid((10,), (-10.0,), (1.0,), (1,))

# expand CartesianGrid with comparable vertices
grid = CartesianGrid((0.0,), (10.0,), dims = (10,))

ll = centroid(grid[1])
uu = centroid(grid[10])

centroid.(grid)
nelements(grid)

left, right = (2,), (3,)
newdim = size(grid) .+ left .+ right

newoffset = offset(grid) .+ left
grid2 = CartesianGrid(newdim, minimum(grid), spacing(grid), newoffset)
@test issubset(vertices(grid), vertices(grid2))

cg = CartesianGrid((-10.0,), (10.0,), dims = (10,))
vertices(cg)

n = 1000
bs = RegularBsplines(-10, 10, n, 3)
x = 0
@time a = basis(x, bs);
@time b = cs.basis_old(x, bs);
a == b

# @benchmark basis(x, bs)
# @benchmark cs.basis_old(x, bs)

n = 10000
bs = RegularBsplines(-10, 10, n, 3)
x = -10:0.001:10
@time a = basis(x, bs);
@time b = cs.basis_sparse(x, bs);
@time c = cs.basis_old(x, bs);
a == b == c
# Matrix(a)

@benchmark basis(x, bs)
@benchmark cs.basis_sparse(x, bs)
# @benchmark cs.basis_old(x, bs)

a == b == c
SparseArrays.getnzval(c)
6002 / (2001*1000)

c[:, 1:1000]

bla = [1, 5, 10]
bla .+ 1
bla .+ 2

vcat([i:i+3 for i in bla]...)

[1:3]

SparseArrays.sparse(repeat(1:11, inner = 3), c[2], c[1])
# reshape(c, 3, length(x))

x = -10:0.01:10
@benchmark a = basis(x, bs)
@benchmark b = cs.basis2(x, bs)
# a == b
#
sparse(rep)

repeat(1:11, inner = 3)

cs.get_x_index(x[1], range(extendedknots(bs)))
map(x -> cs.get_x_index(x, range(extendedknots(bs))), x)

# vcat(collect.([1:3,1:5])...)

x = [0]
c = cs.basis3(x, bs)

# begin
#     for k = 1:10
#         println(k)
#     end
#     print(k)
# end
