# bsplines over regular knots: struct, methods and functions

struct Rbsplines
    minimum::Number
    maximum::Number
    order::Int
    df::Int
end

# df = n_left + 1 + n_internal
function rknots(b::Rbsplines)
    n_left = b.order - 1
    n_right = b.order
    n_internal = b.df - b.order
    return Rknots(b.minimum, b.maximum, n_internal, n_left, n_right)
end

# internal knots
function iknots(b::Rbsplines)
    n_internal = b.df - b.order
    return Rknots(b.minimum, b.maximum, n_internal, 0, 0)
end

function marks(b::Rbsplines)
    knots = knotrange(rknots(b))[1:b.df]
    return knots .+ 0.5 * b.order * step(knots)
end


function rstep(b::Rbsplines)
   return step(knotrange(rknots(b)))
end

function basis(x::Number, b::Rbsplines)

    knots = rknots(b)
    knotset = knotrange(knots)

    # initial arguments
    basis = zeros(length(knotset))

    # order 1
    x_index = get_x_index(x, knotset)
    basis[x_index] = 1

    # order > 1
    δᵣ = zeros(b.order - 1)
    δₗ = zeros(b.order - 1)

    for j = 1:(b.order - 1) # basis order
        δᵣ[j] = knotset[x_index + j] - x
        δₗ[j] = x - knotset[x_index + 1 - j]
        saved = 0
        for r = 1:j # non-zero knots
            term = basis[x_index-j+r] / (δᵣ[r] + δₗ[j+1-r])
            basis[x_index-j-1+r] = saved + δᵣ[r] * term
            saved = δₗ[j+1-r] * term
        end
        basis[x_index] = saved
    end

    return basis[1:b.df]
end

function basis(x::Vector, b::Rbsplines)

    knots = rknots(b)
    knotset = knotrange(knots)

    # initial arguments
    basis = zeros(length(x), length(knotset))

    for i in 1:length(x)
        # order 1
        j_x = get_x_index(x[i], knotset)
        basis[i, j_x] = 1

        # order > 1
        δᵣ = zeros(b.order - 1)
        δₗ = zeros(b.order - 1)

        for j = 1:(b.order - 1) # basis order
            δᵣ[j] = knotset[j_x + j] - x[i]
            δₗ[j] = x[i] - knotset[j_x + 1 - j]
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

function basis(grid::Grid2, b1::Rbsplines, b2::Rbsplines)
    knots1, knots2 = knotrange(grid)
    basis1 = basis(collect(knots1), b1)
    basis2 = basis(collect(knots2), b2)
    return kron(basis1, basis2)
end

function integral(x::Number, b::Rbsplines)
    b = Rbsplines(b.minimum, b.maximum, b.order + 1, b.df + 1)
    ibasis = basis(x, b)
    ibasis = ibasis[2:b.df]
    return rstep(b) * reverse(cumsum(reverse(ibasis)))
end

function integral(x::Vector, b::Rbsplines)
    b = Rbsplines(b.minimum, b.maximum, b.order + 1, b.df + 1)
    ibasis = basis(x, b)
    ibasis = ibasis[:, 2:b.df]
    return rstep(b) * reverse(cumsum(reverse(ibasis, dims = 2), dims = 2), dims = 2)
end

function interval(grid::Vector, b::Rbsplines)
    bs_int = integral(grid, b)
    return diff(bs_int, dims = 1) ./ diff(grid)
end

# function igmrf_marks(order, minimum, maximum, df)
#
#     n_internal = df - order
#     step = (maximum - minimum) / (n_internal + 1)
#
#     maximum = maximum + step * (order)
#     minimum = minimum - step * (order - 1)
#     knots = minimum:step:maximum
#     knots = knots[1:(end-order-1)] .+ step * order / 2
#     return knots
# end
#
# function igmrf_knots(order, minimum, maximum, df)
#
#     n_internal = df - order
#     step = (maximum - minimum) / (n_internal + 1)
#
#     maximum = maximum + step * (order)
#     minimum = minimum - step * (order - 1)
#     knots = minimum:step:maximum + step
#     return knots[1:(end-order-1)]
# end
#
# function grid_markers(min, max, length)
#     step = (max - min) / (length - 1)
#     return (min:step:max)[2:end] .- step / 2
# end
#
