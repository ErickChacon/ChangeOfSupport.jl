# Basic methods for CartesianGrid

"""
    range(x::CartesianGrid)

Return vector of ranges for `x`. Each element of this vector es the range associated to
each dimension of `x`.
"""
function Base.range(x::CartesianGrid{Dim}) where {Dim}
    lower = coordinates(minimum(x))
    upper = coordinates(maximum(x))
    step = spacing(x)
    [(lower[i]):(step[i]):(upper[i]) for i in 1:Dim]
end

"""
    centroids(x::CartesianGrid)

Return vector of centroids for `x`. Each element of this vector is the centroids
associated to each dimension of `x`.
"""
function centroids(x::CartesianGrid{Dim}) where {Dim}
    gridknots = range(x)
    [knots[1:end-1] .+ 0.5 * step(knots) for knots in gridknots]
end

"""
    centroidsmat(x::CartesianGrid)

Return a matrix of centroids for `x`. Each row of this matrix is the centroid associated
to each element of `x`.
"""
function centroidsmat(x::CartesianGrid{Dim}) where {Dim}
    Matrix(transpose(reduce(hcat, coordinates.(centroid.(x)))))
end

"""
    adjacency(x::CartesianGrid; order = 1, cyclic = false)

Return the adjacency `n×n` matrix (A) of the CartesianGrid `x` for a specified `order`.

`n` is the number of elements of the `x`. The value of A[i,j] is `true` if the element `i`
and `j` are neighbors of the specified `order`, and A[i, j] is `false` otherwise. Use the
argument `cyclic = true` to assume neighbors defined by embedding the CartesianGrid in a
torus.
"""
function adjacency(g::CartesianGrid{Dim}; order::Int = 1, cyclic::Bool = false) where {Dim}
    if Dim == 1
        adjacency1(g; order = order, cyclic = cyclic)
    elseif Dim == 2
        adjacency2(g; order = order, cyclic = cyclic)
    else
        throw(ErrorException("adjacency is not implemented for Dim > 2"))
    end
end

function adjacency1(g::CartesianGrid{1}; order::Int = 1, cyclic::Bool = false)
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

function adjacency2(g::CartesianGrid{2}; order::Int = 1, cyclic::Bool = false)
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
