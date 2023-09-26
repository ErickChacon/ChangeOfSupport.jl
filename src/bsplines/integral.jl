"""
    basis(x::CartesianGrid{1}, bs::RegularBsplines)

Returns an sparse evaluation of the basis splines `bs` at elements of `x`.

Each element of `x` is a `Segment`. We consider the evaluation of `bs` at a `Segment` as
mean measure defined as the integral of the basis splines over it divided by the length of
the `Segment`. The output is a sparse matrix of size `length(x)×bs.df`.
"""
function basis(x::CartesianGrid{1}, bs::RegularBsplines)
    gridknots = range(x)[1]
    # 1) compute the integral of the non-zero basis from -∞ to each vertex
    ibasis, indices = nonzerointegral(gridknots, bs)
    # 2) difference between the integral at ending vertices and starting vertices
    output = sparsebasis(ibasis[2:end, :], indices[2:end], bs.df)
    output -= sparsebasis(ibasis[1:(end-1), :], indices[1:(end-1)], bs.df)
    # 3) add constant fullintegral for splines that are complete at ending vertices
    output += fullintegral(indices[1:(end-1)] .- bs.order .+ 1, indices[2:end] .- bs.order, bs)
    # 4) obtain mean measure by dividing integral of segments over length of segments
    output /= step(gridknots)
end

"""
    integral(x, bs::RegularBsplines)

Returns an sparse evaluation of the integral of the bases splines `bs` from `-∞` to `x`.

The returning object is sparse because at any value the integral of the basis splines are
zero for those sarting in knots higher then `x`. The output is a sparse vector of size
`bs.df` if `x` is a `Number`, and is a sparse matrix of size `length(x)×bs.df` if `x` is a
`Vector` or `AbstractRange`.
"""
function integral(x::Union{AbstractRange,Vector}, bs::RegularBsplines)
    # 1) computes the integral of `bs.order` non-zero basis splines and the indices for
    # the last non-zero basis splines of each element in `x`.
    ibasis, lastindices = nonzerointegral(x, bs)
    # 2) converts the integral to a sparse array associated to the `bs.df` basis splines
    # at `x`.
    output = sparsebasis(ibasis, lastindices, bs.df)
    # 3) adds constant term for the full integrals
    output += fullintegral([1], lastindices .- bs.order, bs)
end

"""
    fullintegral(jstart::Vector, jfinish::Vector, bs::RegularBsplines)

Returns a sparse matrix with the full integral value of a basis spline of `bs` from `jstart`
to `jfinish`.

The full integral is simply the step size of the knots associated to `bs`. This
function is used to complete the evaluation of integrals in the case that it is a constant
known value instead of performing computationally expensive operations.
"""
function fullintegral(jstart::Vector, jfinish::Vector, bs::RegularBsplines)
    n = max(length(jstart), length(jfinish))
    I = vcat(map(fill, 1:n, jfinish .- jstart .+ 1)...)
    J = vcat(broadcast(:, jstart, jfinish)...)
    sparse(I, J, step(bs), n, bs.df + 1)[:, 1:bs.df]
end

"""
    nonzerointegral(x::Union{AbstractRange,Vector}, bs::RegularBsplines)

Returns the integral from `-∞` to `x` of the `bs.order` non-zero basis splines at `x` and
indices for the last non-zero basis splines.

The returning matrix has `length(x)` rows and `bs.order` columns, each row correspond to
the intgral of non-zero basis splines of each element in `x`. The indices of the last
non-zero basis splines are also returned to know which basis splines are non-zero. The
algorithm is based on De Boor (2001).
"""
function nonzerointegral(x::Union{AbstractRange,Vector}, bs::RegularBsplines)
    # define new knots and step
    bi = RegularBsplines(bs.lower, bs.upper, bs.df + 1, bs.order + 1)
    knotstep = step(range(extendedknots(bi)))
    # compute non-zero basis of order `bs.order + 1`
    ibasis, indices = nonzerobasis(x, bi)
    # remove first column
    ibasis = ibasis[:, 2:(bi.order)]
    indices = indices .-1
    # compute integral
    ibasis = knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
    return ibasis, indices
end

# ------------------------------
# Old functions. Needs cleaning.
# ------------------------------


"""
    integral_old(x::Number, b::RegularBsplines)

Return an sparse vector of the integral of the basis functions of `b` up to `x`.
"""
function integral_old(x::Number, b::RegularBsplines)
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.df + 1, b.order + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis)))
end

"""
    integral(x::Vector, b::RegularBsplines)

Return an sparse matrix of the integral of the basis functions of `b` up to `x`.
"""
function integral_old(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    # TODO: this implementation is slow for large vectors
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.df + 1, b.order + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[:, 2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
end

# knots where basis functions start and upper boundary
function startingknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, -1)
end

# function boundaryknots(b::RegularBsplines)
#     RegularKnots(b.lower, b.upper, b.df - b.order, 0, 0)
# end


# # centroid reference for basis function
# function centroids(b::RegularBsplines)
#     h = (b.upper - b.lower) / (b.df - b.order + 1)
#     RegularKnots(b.lower + h * b.order / 2,
#                  b.upper + h * b.order / 2,
#                  b.df - b.order, b.order - 1, -1)
# end
#
# # knots for centroid reference for basis function
# function centroidknots(b::RegularBsplines)
#     h = (b.upper - b.lower) / (b.df - b.order + 1)
#     RegularKnots(b.lower + h * (b.order - 1) / 2,
#                  b.upper + h * (b.order - 1) / 2,
#                  b.df - b.order, b.order - 1, 0)
# end




function basis_old(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    knots = extendedknots(b)
    knotrange = range(knots)

    # initial arguments
    basis = zeros(length(x), length(knotrange))

    # check x inside bsplines domain
    if !all(b.lower ≤ val ≤ b.upper for val in x)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    for i in 1:length(x)
        # order 1
        j_x = get_x_index(x[i], knotrange)
        basis[i, j_x] = 1

        # order > 1
        δᵣ = zeros(b.order - 1)
        δₗ = zeros(b.order - 1)

        for j = 1:(b.order - 1) # basis order
            δᵣ[j] = knotrange[j_x + j] - x[i]
            δₗ[j] = x[i] - knotrange[j_x + 1 - j]
            saved = 0
            for r = 1:j # non-zero knots
                term = basis[i, j_x-j+r] / (δᵣ[r] + δₗ[j+1-r])
                basis[i, j_x-j-1+r] = saved + δᵣ[r] * term
                saved = δₗ[j+1-r] * term
            end
            basis[i, j_x] = saved
        end
    end

    return basis[:, 1:b.df]
end

function basis(x::RectilinearGrid{1}, b::RegularBsplines)
    gridknots = knotset(x)[1]
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ diff(gridknots)
end

function basis(x::Segment{1}, b::RegularBsplines)
    gridknots = map(x -> coordinates(x)[1], extrema(x)) |> extrema |> collect
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ diff(gridknots)
end

function basis(x::GeometrySet{1, Float64, Segment{1, Float64}}, b::RegularBsplines)
    vcat([basis(s, b) for s in x]...)
end

# function cartesiangrid(b::RegularBsplines)
#     n_knots = b.df + b.order + 1
#     step = (b.upper - b.lower) / (b.df - b.order + 1)
#     CartesianGrid((n_knots - 1,), (b.lower,), (step,), (b.order,))
#     # RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, b.order)
# end

