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

# get adjacency matrix of CartesianGrid
function adjacency(g::CartesianGrid{1})
    n = g.dims[1]

    # neighbors to the right (→)
    Ir = 1:(n-1)
    Jr = 2:n
    # neighbors to the left (←) are obtained by exchanging the I and J indices.

    # sparse adjacency matrix
    sparse(vcat(Ir, Jr), vcat(Jr, Ir), true, n, n)
end

function adjacency(g::CartesianGrid{1})
    n = g.dims[1]

    # neighbors to the right (→)
    Ir = 1:(n-1)
    Jr = 2:n
    # neighbors to the left (←) are obtained by exchanging the I and J indices.

    # sparse adjacency matrix
    sparse(vcat(Ir, Jr), vcat(Jr, Ir), true, n, n)
end

function adjacency_cyclic(g::CartesianGrid{1})
    n = g.dims[1]

    # neighbors to the right (→)
    Ir = 1:(n-1)
    Jr = 2:n
    Irc = n
    Jrc = 1
    # neighbors to the left (←) are obtained by exchanging the I and J indices.

    # sparse adjacency matrix
    sparse(vcat(Ir, Irc, Jr, Jrc), vcat(Jr, Jrc, Ir, Irc), true, n, n)
end

# get adjacency matrix of CartesianGrid
function adjacency(g::CartesianGrid{2})
    n1, n2 = g.dims
    n = n1 * n2

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        i + (j - 1) * n1
    end

    # neighbors to the right (→)
    Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
    Jr = [ij_to_k(i + 1, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
    # neighbors to the top (↑)
    It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
    Jt = [ij_to_k(i, j + 1, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
    # neighbors to the left (←) and bottom (↓) are obtained by exchanging the I and J
    # indices.

    # sparse adjacency matrix
    sparse(vcat(Ir, It, Jr, Jt), vcat(Jr, Jt, Ir, It), true, n, n)
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
