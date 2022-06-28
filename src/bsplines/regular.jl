"""
    RegularBsplines(lower, upper, order, df)

Bsplines of specified `order` with `df` degrees of freedom, using regular knots with
`lower` and `upper` bounds.

In order to evaluate the basis functions, we add `order - 1` knots to the left of `lower`
such as there are `order - 1` basis starting at the additional left knots, 1 basis
starting at `lower` and the remainin basis starts at the `ni` internal knots (ni = df -
order). We return, `df` basis functions of the specified `order`.
"""
struct RegularBsplines
    lower::Number
    upper::Number
    df::Int
    order::Int
end

"""
    extendedknots(b::RegularBsplines)

Return an extended RegularKnots object necessary to evaluate the RegularBsplines `b`.
"""
function extendedknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, b.order)
end

"""
    CartesianGrid(b::RegularBsplines)

Return a `CartesianGrid` associated with the basis splines of `b`.

The length of `b` is equals to the degrees of freedom of `b`. The returning
`CartesianGrid` is built such as its centroids are the same as the the centroids of the
basis functions of `b`.
"""
function CartesianGrid(b::RegularBsplines)
    spacing = (b.upper - b.lower) /  (b.df - b.order + 1)
    shift = spacing * (b.order - 1) / 2
    CartesianGrid((b.df,), (b.lower + shift,), (spacing,), (b.order,))
    # nknots = b.df + b.order + 1
    # CartesianGrid((nknots-1,), (b.lower,), (spacing,), (b.order,))
end

"""
    get_x_index(x::Number, knots::AbstractRange)

Get index i such as `knots[i] ≤ x < knots[i+1]`. We use it to identify basis splines of
certain order that are non-zero at x.
"""
function get_x_index(x::Number, knots::AbstractRange)

    # check if x is inside the range
    if x < knots[1]
        return 0
    elseif x >= knots[end]
        return length(knots)
    end

    # compute index if x is inside the range
    i = 1 + floor(Int, (x - knots[1]) / step(knots))

    # simple fix in case condition knots[i] ≤ x < knots[i+1] is not hold
    if x >= knots[i + 1]
        i = i + 1
    end
    if x < knots[i]
        i = i - 1
    end

    return i
end


"""
    basis(x::Number, b::RegularBsplines)

Return an sparse vector of the evaluation of the basis functions of `b` at `x`.

The dimension of the returning sparse vector is `b.df`. The algorithm is based on De Boor
(2001, page 110).
"""
function basis(x::Number, b::RegularBsplines)
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis
    basis = spzeros(length(knotrange))

    # check x inside acceptable range
    if !(b.lower ≤ x ≤ b.upper)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # order = 1
    j = get_x_index(x, knotrange)
    basis[j] = 1

    # order > 1
    δᵣ = zeros(b.order - 1)
    δₗ = zeros(b.order - 1)

    # uses basis of order k to compute basis of order k+1
    for k = 1:(b.order - 1)
        # create left and right differences
        δᵣ[k] = knotrange[j + k] - x
        δₗ[k] = x - knotrange[j + 1 - k]
        leftterm = 0
        # evaluate the first k non-zero basis of order k+1 (s = 1:k)
        for s = 1:k
            zs = basis[j-k+s] / (δᵣ[s] + δₗ[k-s+1])
            basis[j-k-1+s] = leftterm + δᵣ[s] * zs
            # leftterm for s' = s+1, then δ(k-s'+2) = δ(k-s+1) and z(s'-1) = z(s).
            leftterm = δₗ[k-s+1] * zs
        end
        # evaluate the k+1 basis of order k+1 (s = k+1)
        basis[j] = leftterm
    end
    # return only the first df basis functions, these are non-zero in the desired range
    basis[1:b.df]
end

"""
    basis(x::AbstractVector, b::RegularBsplines)

Return an sparse matrix of the evaluation of Bsplines `b` at `x`.

The sparse matrix has `length(x)` rows and `b.df` columns. The algorithm is based on De
Boor (2001, page 110).
"""
function basis(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis and indices to create sparse arrays (n*K size)
    basis = zeros(length(x) * b.order) # basis values
    I = repeat(1:length(x), inner = b.order) # index of observation
    J = zeros(Int, length(x) * b.order) # index of basis function

    # check x inside bsplines domain
    if !all(b.lower ≤ val ≤ b.upper for val in x)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    for i in 1:length(x)
        # order = 1
        ibase = (i-1) * b.order
        j_x = get_x_index(x[i], knotrange)
        basis[ibase + 1] = 1

        # order > 1
        δᵣ = zeros(b.order - 1)
        δₗ = zeros(b.order - 1)

        # uses basis of order k to compute basis of order k+1
        for k = 1:(b.order - 1)
            # create left and right differences
            δᵣ[k] = knotrange[j_x + k] - x[i]
            δₗ[k] = x[i] - knotrange[j_x + 1 - k]
            leftterm = 0
            # evaluate the first k non-zero basis of order k+1 (s = 1:k)
            for s = 1:k
                zs = basis[ibase + s] / (δᵣ[s] + δₗ[k+1-s])
                basis[ibase + s] = leftterm + δᵣ[s] * zs
                # leftterm for s' = s+1, then δ(k-s'+2) = δ(k-s+1) and z(s'-1) = z(s).
                leftterm = δₗ[k+1-s] * zs
            end
            # evaluate the k+1 basis of order k+1 (s = k+1)
            basis[ibase + k + 1] = leftterm
        end

        # J indices for non-zero basis
        J[(ibase+1):(ibase+b.order)] = (j_x-b.order+1):j_x
    end

    # return only the first df basis functions
    sparse(I, J, basis, length(x), length(knotrange))[:, 1:b.df]
end

"""
    integral(x::Number, b::RegularBsplines)

Return an sparse vector of the integral of the basis functions of `b` up to `x`.
"""
function integral(x::Number, b::RegularBsplines)
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
function integral(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.df + 1, b.order + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[:, 2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
end




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

