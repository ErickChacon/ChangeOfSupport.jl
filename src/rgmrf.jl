"""
    RGMRF(grid, order, δ, κ)

Regular Gaussian Markov random field with n-th order, δ as and κ as precision.
"""
struct RGMRF{N,T}
    grid::CartesianGrid{N,T}
    order::Integer
    δ::Number
    κ::Number
end

function structure(X::RGMRF{N,T}) where {N,T}
    A = adjacency(X.grid)
    if X.order == 1
        neighs = vec(sum(A, dims = 2))
        Diagonal(neighs .+ X.δ) - A
    else
        throw(ErrorException("not implemented"))
    end
end

function Base.rand(X::RGMRF{N,T}) where {N,T}
    n = nelements(X.grid)
    R = structure(X)
    x = cholesky(R).UP \ randn(n)
    x / sqrt(X.κ)
end

