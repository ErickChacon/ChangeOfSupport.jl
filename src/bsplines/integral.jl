"""
    integral(x::Number, b::RegularBsplines)

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
    # TODO: this implimentation is slow for large vectors
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.df + 1, b.order + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[:, 2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
end

function _integral(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    # define new knots and step
    bi = RegularBsplines(b.lower, b.upper, b.df + 1, b.order + 1)
    knotstep = step(range(extendedknots(bi)))
    # compute integral
    ibasis, indices = nonzerobasis(x, bi)
    # remove first column
    ibasis = ibasis[:, 2:(bi.order)]
    indices = indices .-1
    # compute integral
    ibasis = knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
    return ibasis, indices
end

function integral(x::Union{AbstractRange,Vector}, bs::RegularBsplines)
    # 1) we compute only the `bs.order` non-zero basis splines and the indices for the
    # last nons-zero basis spline of each element in `x`.
    ibasis, lastindices = _integral(x, bs)
    # 2) we convert the basis to a sparse design matrix associated to the `bs.df` basis
    # splines.
    sparsebasis(ibasis, lastindices, bs.df)
end

# function basis(x::CartesianGrid{1}, b::RegularBsplines)
#     gridknots = range(x)[1]
#     ibasis = integral(gridknots, b)
#     diff(ibasis, dims = 1) ./ step(gridknots)
# end


# function boundaryknots(b::RegularBsplines)
#     RegularKnots(b.lower, b.upper, b.df - b.order, 0, 0)
# end

# knots where basis functions start and upper boundary
function startingknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, -1)
end

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


# # algorithm based on boor 2001, page 110
# function basis_old(x::Number, b::RegularBsplines)
#     knots = extendedknots(b)
#     knotrange = range(knots)
#
#     # initial arguments
#     basis = zeros(length(knotrange))
#
#     # check x inside acceptable range
#     if !(b.lower ≤ x ≤ b.upper)
#         throw(DomainError(x, "The values of x should lie in the range of b."))
#     end
#
#     # order = 1
#     x_index = get_x_index(x, knotrange)
#     basis[x_index] = 1
#
#     # order > 1
#     δᵣ = zeros(b.order - 1)
#     δₗ = zeros(b.order - 1)
#
#     for j = 1:(b.order - 1) # basis order
#         δᵣ[j] = knotrange[x_index + j] - x
#         δₗ[j] = x - knotrange[x_index + 1 - j]
#         saved = 0
#         for r = 1:j # non-zero knots
#             term = basis[x_index-j+r] / (δᵣ[r] + δₗ[j+1-r])
#             basis[x_index-j-1+r] = saved + δᵣ[r] * term
#             saved = δₗ[j+1-r] * term
#         end
#         basis[x_index] = saved
#     end
#
#     basis[1:b.df]
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

# function basis_sparse(x::Union{AbstractRange,Vector}, b::RegularBsplines)
#     knots = extendedknots(b)
#     knotrange = range(knots)
#
#     # initial arguments
#     basis = spzeros(length(x), length(knotrange))
#
#     # check x inside bsplines domain
#     if !all(b.lower ≤ val ≤ b.upper for val in x)
#         throw(DomainError(x, "The values of x should lie in the range of b."))
#     end
#
#     for i in 1:length(x)
#         # order = 1
#         j_x = get_x_index(x[i], knotrange)
#         basis[i, j_x] = 1
#
#         # order > 1
#         δᵣ = zeros(b.order - 1)
#         δₗ = zeros(b.order - 1)
#
#         for k = 1:(b.order - 1) # basis order
#             δᵣ[k] = knotrange[j_x + k] - x[i]
#             δₗ[k] = x[i] - knotrange[j_x + 1 - k]
#             leftterm = 0
#             for s = 1:k # non-zero knots
#                 zs = basis[i, j_x-k+s] / (δᵣ[s] + δₗ[k+1-s])
#                 basis[i, j_x-k-1+s] = leftterm + δᵣ[s] * zs
#                 leftterm = δₗ[k+1-s] * zs
#             end
#             basis[i, j_x] = leftterm
#         end
#     end
#
#     return basis[:, 1:b.df]
# end


function basis(x::CartesianGrid{1}, b::RegularBsplines)
    gridknots = range(x)[1]
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ step(gridknots)
end

function basis(x::RectilinearGrid{1}, b::RegularBsplines)
    gridknots = knotset(x)[1]
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ diff(gridknots)
end

# function basis(x::Segment{1}, b::RegularBsplines)
#     # gridknots = knotset(x)[1]
#     # ibasis = integral(gridknots, b)
#     # diff(ibasis, dims = 1) ./ diff(gridknots)
# end

# function cartesiangrid(b::RegularBsplines)
#     n_knots = b.df + b.order + 1
#     step = (b.upper - b.lower) / (b.df - b.order + 1)
#     CartesianGrid((n_knots - 1,), (b.lower,), (step,), (b.order,))
#     # RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, b.order)
# end

