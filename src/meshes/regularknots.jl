"""
    RegularKnots(lower, upper, ni, nl, nr)

Defines a RegularKnots struct having `ni` knots between `lower` and `upper`, and with
additional `nl` knots at the left and `nr` knots at the right. In total, it has `ni + nl +
nr + 2` knots because it includes `lower` and `upper` as knots.
"""
struct RegularKnots
    lower::Float64
    upper::Float64
    ni::Int
    nl::Int
    nr::Int
end

"""
Number of knots in a `RegularKnots` object.
"""
Base.length(x::RegularKnots) = x.ni + x.nl + x.nr + 2

"""
Converts a `RegularKnots` object `x` to a `range` using `x.lower` as reference value. It
allows to create ranges with comparable knots because it uses the same reference value.
"""
function Base.range(x::RegularKnots)
    n_knots = length(x)
    step = (x.upper - x.lower) / (x.ni + 1)
    knotrange = StepRangeLen(
        Base.TwicePrecision{Float64}(x.lower),
        Base.TwicePrecision{Float64}(step),
        n_knots, x.nl + 1)
    return knotrange
end

# RegularKnots are important to compute bsplines because it allows to build knots based on
# a reference value and later extend it to obtain numerically comparable knots. In case we
# do not use a unique reference value, the knots are not always comparable in practice due
# to limited precision for real numbers. In case of bsplines, the recursive relationship,
# and also the recursive integral relationship need to work with a set of knots that are
# reduced or extended depending of the of the `order` of the bsplines.
