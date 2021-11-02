function Base.range(x::CartesianGrid{Dim}) where {Dim}
    lower = coordinates(minimum(x))
    upper = coordinates(maximum(x))
    step = spacing(x)
    [(lower[i]):(step[i]):(upper[i]) for i in 1:Dim]
end

function centroids(x::CartesianGrid{Dim}) where {Dim}
    gridknots = range(x)
    [knots[1:end-1] .+ 0.5 * step(knots) for knots in gridknots]
end

struct IrregularGrid
    minimum::Number
    maximum::Number
    iknots::Vector
end

function Meshes.vertices(x::IrregularGrid)
    [x.minimum, sort(x.iknots)..., x.maximum]
end

function centroids(x::IrregularGrid)
    knots = vertices(x)
    knots[1:end-1] + diff(knots) / 2
end

# function marks(x::Igrid)
#     knots = knotset(x)
#     return knots[1:end-1] + diff(knots) / 2
# end


# function centroids(x::CartesianGrid{1})
#     lower = coordinates(minimum(x))[1]
#     upper = coordinates(maximum(x))[1]
#     step = spacing(x)[1]
#     lower:step:upper
# end


# struct Grid
#     minimum::Number
#     maximum::Number
#     df::Int
# end
#
# function knotrange(x::Grid)
#     step = (x.maximum - x.minimum) / (x.df)
#     knotrange = StepRangeLen(
#         Base.TwicePrecision{Float64}(x.minimum),
#         Base.TwicePrecision{Float64}(step),
#         x.df + 1)
#     return knotrange
# end
#
# function marks(x::Grid)
#     knots = knotrange(x)
#     return knots[1:end-1] .+ 0.5 * step(knots)
# end
#
#
# struct Grid2
#     minimum1::Number
#     maximum1::Number
#     df1::Int
#     minimum2::Number
#     maximum2::Number
#     df2::Int
# end
#
# function knotrange(x::Grid2)
#     g1 = Grid(x.maximum1, x.minimum1, x.df1)
#     g2 = Grid(x.maximum2, x.minimum2, x.df2)
#     return knotrange(g1), knotrange(g2)
# end
#
# function knotset(x::Grid2)
#     knots1, knots2 = knotrange(x)
#     s1 = repeat(knots1, inner = length(knots2))
#     s2 = repeat(knots2, outer = length(knots1))
#     return [s1 s2]
# end
#
# # knotset(Grid2(-10, 10, 3, -5, 5, 4))
