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

"""
Return the adjacency `n×n` matrix (A) of a `CartesianGrid` where `n` is the number of
elements of the `CartesianGrid`. The value of A[i,j] is 1 if the element `i` and `j` are
neighbors of the specified `order`, and A[i, j] is 0 otherwise. Use the argument `cyclic =
true` to assume neighbors defined by embedding the CartesianGrid in a torus.
"""
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
