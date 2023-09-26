"""
    knotset(x::RectilinearGrid{Dim, T})

Return a vector of size `Dim` of knot sets. Each element of this vector is the knot set of
each dimension of `x`.
"""
knotset(d::Meshes.RectilinearGrid) = d.xyz

"""
    centroids(x::RectilinearGrid{Dim,T})

Return a vector of size `Dim` of centroids. Each element of this vector is the centroids of
each dimension of `x`.
"""
function centroids(x::Meshes.RectilinearGrid{Dim,T}) where {Dim,T}
    knots = knotset(x)
    [knots[i][1:end-1] + diff(knots[i]) / 2 for i in 1:Dim]
end

# Get vector of marginal centroids for each dimension of a `RectilinearGrid` object.

function centroids(x::GeometrySet{1, Float64, Segment{1, Float64}})
    (map(s -> coordinates(centroid(s))[1], x),)
end

function knotset(x::GeometrySet{1, Float64, Segment{1, Float64}})
    (vcat(map(s -> coordinates(minimum(s))[1], x), map(s -> coordinates(maximum(s))[1], x)),)
end
