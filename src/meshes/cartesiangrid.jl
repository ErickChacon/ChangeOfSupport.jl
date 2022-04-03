## Basic methods for CartesianGrid

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

## Adjacency sparse matrix of CartesianGrid (temporal, spatial)

function adjacency(g::CartesianGrid{1}; order = 1, cyclic = false)
    n = nelements(g)

    # neighbors to the right (→)
    if cyclic
        Ir = 1:n
        Jr = mod1.(Ir .+ order, n)
    else
        Ir = 1:(n - order)
        Jr = (1 + order):n
    end

    # sparse adjacency matrix
    A = sparse(Ir, Jr, true, n, n)
    A .| A'
end

function adjacency(g::CartesianGrid{2}; order = 1, cyclic = false)
    n1, n2 = size(g)
    n = nelements(g)

    # function to convert i,j cel to k-index (CartesianGrid)
    function ij_to_k(i, j, n1, n2)
        (j - 1) * n1 + i
    end

    if cyclic
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:n2]
        Jr = [ij_to_k(mod1(i + order, n1), j, n1, n2) for i = 1:n1 for j = 1:n2]
        # neighbors to the top (↑)
        It = Ir
        Jt = [ij_to_k(i, mod1(j + order, n2), n1, n2) for i = 1:n1 for j = 1:n2]
    else
        # neighbors to the right (→)
        Ir = [ij_to_k(i, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        Jr = [ij_to_k(i + order, j, n1, n2) for i = 1:(n1-order) for j = 1:n2]
        # neighbors to the top (↑)
        It = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
        Jt = [ij_to_k(i, j + order, n1, n2) for i = 1:n1 for j = 1:(n2-order)]
    end

    # sparse adjacency matrix
    A = sparse(vcat(Ir, It), vcat(Jr, Jt), true, n, n)
    A .| A'
end

## Differences and structure for GMRF

function difference(g::CartesianGrid{1}; order = 1, cyclic = false)
    n = nelements(g)

    if order == 1 && !cyclic
        return spdiagm(n-1, n, 0 => fill(-1, n-1), 1 => fill(1, n-1))
    end

    if order == 1 && cyclic
        return spdiagm(0 => fill(-1, n), 1 => fill(1, n-1), -(n-1) => fill(1, 1))
    end

    if order == 2 && !cyclic
        return spdiagm(n-2, n,
                       0 => fill(1, n-2), 1 => fill(-2, n-2), 2 => fill(1, n-2))
    end

    if order == 2 && cyclic
        return spdiagm(0 => fill(1, n), 1 => fill(-2, n-1), 2 => fill(1, n-2),
                       -(n-2) => fill(1, 2), -(n-1) => fill(-2, 1))
    end

    throw(ErrorException("not implemented"))
end

function structure(g::CartesianGrid{1}; δ = 0, order = 1, cyclic = false)
    D = difference(g, order = order, cyclic = cyclic)
    D'D + δ * I
end
