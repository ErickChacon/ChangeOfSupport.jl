"""
    RegularBsplines(lower, upper, order, df)

Basis splines of specified `order` with `df` degrees of freedom, using regular knots with
`lower` and `upper` bounds.

The `order` is the number of parameters required to represent the polinomial function.
For example, a quadratic polinomial has `order` 3. The degrees of freedom `df` is the
number of desired basis splines.

In order to evaluate the `df` basis functions using the recursive relationship of basis
splines. We need to create `order - 1` knots to the left of `lower` which are the starting
points of `order - 1` basis functions. One basis function starts at `lower` and the
remaining functions start at the `df - order` internal knots.
"""
struct RegularBsplines
    lower::Number
    upper::Number
    df::Int
    order::Int
end

Base.minimum(bs::RegularBsplines) = bs.lower
Base.maximum(bs::RegularBsplines) = bs.upper
Base.extrema(bs::RegularBsplines) = minimum(bs), maximum(bs)
Base.length(bs::RegularBsplines) = bs.df
order(bs::RegularBsplines) = bs.order

# ----------------------------------------------
# Auxiliary functions to evaluate basis splines.
# ----------------------------------------------

"""
    extendedknots(bs::RegularBsplines)

Return an extended `RegularKnots` object necessary to evaluate the RegularBsplines `bs`
using the recurrent relationship.

Creates a sequence of `knots` where there are `bs.order - 1` knots to the left of
`bs.lower`, `bs.order` knots to the right of `bs.upper` and there are `bs.df - bs.order`
internal knots between `bs.left` and `bs.right`.
"""
function extendedknots(bs::RegularBsplines)
    RegularKnots(bs.lower, bs.upper, bs.df - bs.order, bs.order - 1, bs.order)
end

step(bs::RegularBsplines) = step(range(extendedknots(bs)))

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

