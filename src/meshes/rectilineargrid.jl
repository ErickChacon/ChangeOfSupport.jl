"""
    RectilinearGrid(lower, upper, iknots)

Define a `RectilinearGrid` with `iknots` knots between `lower` and `upper`.
"""
struct RectilinearGrid{Dim,T}
    lower::NTuple{Dim,T}
    upper::NTuple{Dim,T}
    iknots::NTuple{Dim, Vector{T}}

    function RectilinearGrid{Dim,T}(lower, upper, iknots) where {Dim,T}
        iknots = map(unique ∘ sort, iknots)
        new{Dim,T}(lower,upper,iknots)
    end
end

function RectilinearGrid(lower::NTuple{Dim,T}, upper::NTuple{Dim,T},
                         iknots::NTuple{Dim, Vector{T}}) where {Dim,T}
    RectilinearGrid{Dim,T}(lower, upper, iknots)
end

"""
    knotset(x::RectilinearGrid{Dim, T})

Return a vector of size `Dim` of knot sets. Each element of this vector is the knot set of
each dimension of `x`.
"""
function knotset(x::RectilinearGrid{Dim,T}) where {Dim,T}
    [[x.lower[i], x.iknots[i]..., x.upper[i]] for i in 1:Dim]
end

"""
    centroids(x::RectilinearGrid{Dim,T})

Return a vector of size `Dim` of centroids. Each element of this vector is the centroids of
each dimension of `x`.
"""
function centroids(x::RectilinearGrid{Dim,T}) where {Dim,T}
    knots = knotset(x)
    [knots[i][1:end-1] + diff(knots[i]) / 2 for i in 1:Dim]
end
# Get vector of marginal centroids for each dimension of a `RectilinearGrid` object.
