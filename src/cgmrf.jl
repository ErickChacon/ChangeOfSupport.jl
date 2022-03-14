"""
    CGMRF(grid, order, δ, κ)

Circulant Gaussian Markov random fields.
"""
struct CGMRF{N,T}
    grid::CartesianGrid{N,T}
    order::Integer
    δ::Number
    κ::Number
end

function structure_base(g::CartesianGrid{1}; δ = 0, order = 1)
    n = g.dims[1]
    base = zeros(n)
    if order == 1
        base[[1, 2, end]] = [2 + δ, -1, -1]
    elseif order == 2
        base[[1, 2, 3, end-1, end]] = [6 + δ, -4, 1, 1, -4]
    else
        throw(ErrorException("not implemented"))
    end
    base
end

function structure_base(g::CartesianGrid{2}; δ = 0, order = 1)
    (n1, n2) = g.dims
    n = n1 * n2
    base = zeros(n2, n1)
    if order == 1
        base[1, 1] = 4 + δ
        base[1, 2] = -1
        base[2, 1] = -1
        base[1, end] = -1
        base[end, 1] = -1
    elseif order == 2
        base[1, 1] = 20 + δ
        base[1, 2] = -8
        base[2, 1] = -8
        base[1, end] = -8
        base[end, 1] = -8
        base[2, 2] = 2
        base[2, end] = 2
        base[end, 2] = 2
        base[end, end] = 2
        base[1, 3] = 1
        base[3, 1] = 1
        base[1, end-1] = 1
        base[end-1, 1] = 1
    else
        throw(ErrorException("not implemented"))
    end
    base
end

function structure(X::CGMRF{1,T}) where {T}
    A = adjacency(X.grid)
    if X.order == 1
        (2 + X.δ)I - A
    end
end


function Base.rand(X::CGMRF{N,T}; border = 0) where {N,T}
    dims = X.grid.dims
    n = prod(dims)

    # define base
    base = structure_base(X.grid; δ = X.δ, order = X.order)

    # compute eigenvalues
    λ = fft(base)

    # simulate
    z = randn(reverse(dims)...) + randn(reverse(dims)...) * im
    x = real(fft(λ .^ (-1/2) .* z) / sqrt(n))
    x = x / sqrt(X.κ)

    # edge effects
    if length(dims) == 1
        x[border + 1:n - border]
    elseif length(dims) == 2
        x[border + 1:dims[2] - border, border + 1:dims[1] - border]
    end
end
