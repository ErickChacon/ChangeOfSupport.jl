"""
    NRegularBsplines(lower, upper, df, order)

N-dimensional tensor product basis splines of specified `order` with `df` degrees of
freedom, using regular knots with `lower` and `upper` bounds.
"""
struct NRegularBsplines{N,T}
    lower::NTuple{N,T}
    upper::NTuple{N,T}
    df::NTuple{N, Int}
    order::Int
end

function NRegularBsplines(grid::CartesianGrid{N, T}, order::Int) where {N, T}
    lower = Tuple(x for x in coordinates(minimum(grid)))
    upper = Tuple(x for x in coordinates(maximum(grid)))
    NRegularBsplines(lower, upper, size(grid), order)
end

# ----------------------------------------------
# Sparse evaluation of NRegularBsplines.
# ----------------------------------------------

function basis(x::NTuple{N, T}, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    bb = [basis(x[i], b[i]) for i in 1:N]
    N == 1 ? bb[1] : kron(reverse(bb)...)
end

# this need to be faster
function basis(x::Array{T, N}, b::NRegularBsplines{N, T}) where {N, T}
    B = zeros(size(x)[1], prod(b.df))
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    bb = [basis(x[:, i], b[i]) for i in N:-1:1]
    @time N == 1 ? bb[1] : kron.(eachrow.(bb)...) |> x -> Matrix(transpose(reduce(hcat, x)))
end

# this need to be faster
function basis2(x::Array{T, N}, b::NRegularBsplines{N, T}) where {N, T}
    B = zeros(size(x)[1], prod(b.df))
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    bb = [basis(x[:, i], b[i]) for i in N:-1:1]
    # @time N == 1 ? bb[1] : kron.(eachrow.(bb)...) |> x -> Matrix(transpose(reduce(hcat, x)))
end


function basis(x::CartesianGrid{N}, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    gridknots = range(x)
    ibasis = [integral(gridknots[i], b[i]) for i in 1:N]
    ibasis = [diff(ibasis[i], dims = 1) ./ step(gridknots[i]) for i in 1:N]
    N == 1 ? ibasis[1] : kron(reverse(ibasis)...)
end

function basis(x::RectilinearGrid{N}, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    gridknots = knotset(x)
    ibasis = [integral(gridknots[i], b[i]) for i in 1:N]
    ibasis = [diff(ibasis[i], dims = 1) ./ diff(gridknots[i]) for i in 1:N]
    N == 1 ? ibasis[1] : kron(reverse(ibasis)...)
end

function basis(x::Quadrangle{N}, b::NRegularBsplines{N, T}) where {N, T}
    b = [RegularBsplines(b.lower[i], b.upper[i], b.df[i], b.order) for i in 1:N]
    segs = [Segment(Point(coordinates(extrema(x)[1])[i],), Point(coordinates(extrema(x)[2])[i],)) for i = 1:2]
    ibasis = [basis(segs[i], b[i]) for i in 1:N]
    N == 1 ? ibasis[1] : kron(reverse(ibasis)...)
end

function basis(x::GeometrySet{2, Float64, Quadrangle{2, Float64}}, b::NRegularBsplines{N, T}) where {N, T}
    vcat([basis(s, b) for s in x]...)
end

