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

# Converts a RegularKnots object to a range using `x.lower` as reference value. It allows
# to create ranges with comparable knots because it uses the same reference value.
function Base.range(x::RegularKnots)
    n_knots = x.ni + x.nl + x.nr + 2
    step = (x.upper - x.lower) / (x.ni + 1)
    knotrange = StepRangeLen(
        Base.TwicePrecision{Float64}(x.lower),
        Base.TwicePrecision{Float64}(step),
        n_knots, x.nl + 1)
    return knotrange
end

# Get index i such as knots[i] ≤ x < knots[i+1]. We use it to identify basis splines of
# order k that are non-zero at x.
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

# RegularKnots are important to compute bsplines because it allows to build knots based on
# a reference value and later extend it to obtain numerically comparable knots. In case we
# do not use a unique reference value, the knots are not always comparable in practice due
# to limited precision for real numbers. In case of bsplines, the recursive relationship,
# and also the recursive integral relationship need to work with a set of knots that are
# reduced or extended depending of the of the `order` of the bsplines.
