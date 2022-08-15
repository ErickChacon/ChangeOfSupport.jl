"""
    RegularBsplines(lower, upper, order, df)

Basis splines of specified `order` with `df` degrees of freedom, using regular knots with
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

# ----------------------------------------------
# Auxiliary functions to evaluate basis splines.
# ----------------------------------------------

"""
    extendedknots(bs::RegularBsplines)

Return an extended RegularKnots object necessary to evaluate the RegularBsplines `bs`.
"""
function extendedknots(bs::RegularBsplines)
    RegularKnots(bs.lower, bs.upper, bs.df - bs.order, bs.order - 1, bs.order)
end

"""
    get_x_index(x::Number, knots::AbstractRange)

Get index `s` such as `knots[s] ≤ x < knots[s+1]`. We use it to identify basis splines of
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

# ---------------------------------------------------------------------
# Auxiliary function to associate a basis splines with a CartesianGrid.
# ---------------------------------------------------------------------

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

