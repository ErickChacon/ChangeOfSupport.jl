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

function centroidsmat(x::CartesianGrid{Dim}) where {Dim}
    Matrix(transpose(reduce(hcat, coordinates.(centroid.(x)))))
end


struct IrregularGrid{Dim,T}
    minimum::NTuple{Dim,T}
    maximum::NTuple{Dim,T}
    iknots::NTuple{Dim, Vector{T}}

    function IrregularGrid{Dim,T}(minimum, maximum, iknots) where {Dim,T}
        iknots = map(unique ∘ sort, iknots)
        new{Dim,T}(minimum,maximum,iknots)
    end
end

function IrregularGrid(minimum::NTuple{Dim,T}, maximum::NTuple{Dim,T},
                         iknots::NTuple{Dim, Vector{T}}) where {Dim,T}
    IrregularGrid{Dim,T}(minimum, maximum, iknots)
end

function knotset(x::IrregularGrid{Dim,T}) where {Dim,T}
    [[x.minimum[i], x.iknots[i]..., x.maximum[i]] for i in 1:Dim]
end

function centroids(x::IrregularGrid{Dim,T}) where {Dim,T}
    knots = knotset(x)
    [knots[i][1:end-1] + diff(knots[i]) / 2 for i in 1:Dim]
end
