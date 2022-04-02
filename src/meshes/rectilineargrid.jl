"""
    RectilinearGrid(lower, upper, iknots)

Defines a RectilinearGrid struct with `iknots` knots between `lower` and `upper`.
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

function knotset(x::RectilinearGrid{Dim,T}) where {Dim,T}
    [[x.lower[i], x.iknots[i]..., x.upper[i]] for i in 1:Dim]
end

function centroids(x::RectilinearGrid{Dim,T}) where {Dim,T}
    knots = knotset(x)
    [knots[i][1:end-1] + diff(knots[i]) / 2 for i in 1:Dim]
end
