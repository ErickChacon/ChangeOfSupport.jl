"""
    RegularKnots(lower, upper, ni, nl, nr)

Defines a RegularKnots struct having `ni` knots between `lower` and `upper`, and with
additional `nl` knots at the left and `nr` knot at the right.
"""
struct RegularKnots
    lower::Float64
    upper::Float64
    ni::Int
    nl::Int
    nr::Int
end

# build knots as a StepRangeLen type
function Base.range(x::RegularKnots)
    n_knots = x.ni + x.nl + x.nr + 2
    step = (x.upper - x.lower) / (x.ni + 1)
    knotrange = StepRangeLen(
        Base.TwicePrecision{Float64}(x.lower),
        Base.TwicePrecision{Float64}(step),
        n_knots, x.nl + 1)
    return knotrange
end

# get index i such as knots[i] ≤ x < knots[i+1]. This is mainly used to compute bsplines.
function get_x_index(x::Number, knots::AbstractRange)

    # check if x is inside the range
    if x < knots[1]
        return 0
    elseif x >= knots[end]
        return length(knots)
    end

    # compute index is x is inside the range
    i = 1 + floor(Int, (x - knots[1]) / step(knots))

    # simple fix in case it fails because of precision
    if x >= knots[i + 1]
        i = i + 1
    end
    if x < knots[i]
        i = i - 1
    end

    return i
end
