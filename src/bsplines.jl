"""
    RegularBsplines(lower, upper, order, df)
Bsplines defined over a regular grid with `lower` and `upper` bounds. In total we create
`df` basis functions of a specific `order`. In order to evaluate the basis functions, we
add `order - 1` knots to the left of lower such as there are `order - 1` basis starting at
the additional left knots, 1 basis starting at `lower` and the remainin basis starts at
the `ni` internal knots (ni = df - order).
"""
struct RegularBsplines
    lower::Number
    upper::Number
    order::Int
    df::Int
end

function extendedknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, b.order - 1, b.order)
end

function boundaryknots(b::RegularBsplines)
    RegularKnots(b.lower, b.upper, b.df - b.order, 0, 0)
end

function basis(x::Number, b::RegularBsplines)
    knots = extendedknots(b)
    knotrange = range(knots)

    # initial arguments
    basis = zeros(length(knotrange))

    # check x inside acctable range
    if !(b.lower ≤ x ≤ b.upper)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # order 1
    x_index = get_x_index(x, knotrange)
    basis[x_index] = 1

    # order > 1
    δᵣ = zeros(b.order - 1)
    δₗ = zeros(b.order - 1)

    for j = 1:(b.order - 1) # basis order
        δᵣ[j] = knotrange[x_index + j] - x
        δₗ[j] = x - knotrange[x_index + 1 - j]
        saved = 0
        for r = 1:j # non-zero knots
            term = basis[x_index-j+r] / (δᵣ[r] + δₗ[j+1-r])
            basis[x_index-j-1+r] = saved + δᵣ[r] * term
            saved = δₗ[j+1-r] * term
        end
        basis[x_index] = saved
    end

    basis[1:b.df]
end

function basis(x::Union{AbstractRange,Vector}, b::RegularBsplines)

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

# function basis(grid::Grid2, b1::RegularBsplines, b2::RegularBsplines)
#     knots1, knots2 = range(grid)
#     basis1 = basis(collect(knots1), b1)
#     basis2 = basis(collect(knots2), b2)
#     return kron(basis1, basis2)
# end

function integral(x::Number, b::RegularBsplines)
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.order + 1, b.df + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis)))
end

function integral(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    # define new knots and step
    b = RegularBsplines(b.lower, b.upper, b.order + 1, b.df + 1)
    knotstep = step(range(extendedknots(b)))
    # compute integral
    ibasis = basis(x, b)[:, 2:b.df]
    knotstep * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
end

function integral(x::CartesianGrid{1}, b::RegularBsplines)
    gridknots = range(x)
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ step(gridknots)
end

function integral(x::IrregularGrid, b::RegularBsplines)
    gridknots = vertices(x)
    ibasis = integral(gridknots, b)
    diff(ibasis, dims = 1) ./ diff(gridknots)
end






# function interval(grid::Vector, b::RegularBsplines)
#     bs_int = integral(grid, b)
#     return diff(bs_int, dims = 1) ./ diff(grid)
# end
#
# # function igmrf_marks(order, lower, upper, df)
# #
# #     n_internal = df - order
# #     step = (upper - lower) / (n_internal + 1)
# #
# #     upper = upper + step * (order)
# #     lower = lower - step * (order - 1)
# #     knots = lower:step:upper
# #     knots = knots[1:(end-order-1)] .+ step * order / 2
# #     return knots
# # end
# #
# # function igmrf_knots(order, lower, upper, df)
# #
# #     n_internal = df - order
# #     step = (upper - lower) / (n_internal + 1)
# #
# #     upper = upper + step * (order)
# #     lower = lower - step * (order - 1)
# #     knots = lower:step:upper + step
# #     return knots[1:(end-order-1)]
# # end
# #
# # function grid_markers(min, max, length)
# #     step = (max - min) / (length - 1)
# #     return (min:step:max)[2:end] .- step / 2
# # end
# #
#
# function marks(b::RegularBsplines)
#     knots = range(rknots(b))[1:b.df]
#     return knots .+ 0.5 * b.order * step(knots)
# end
#
#
# # function rstep(b::RegularBsplines)
# #    return step(range(rknots(b)))
# # end
#
