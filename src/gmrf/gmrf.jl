"""
    GMRF(S, κ)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: ContinuousMultivariateDistribution
    S::AbstractMatrix
    κ::Real
end

Base.length(d::GMRF) = size(d.S, 1)

function Distributions._rand!(rng::AbstractRNG, d::GMRF, x::AbstractVector{T}) where T<:Real
    randn!(rng, x)
    copyto!(x, cholesky(d.S).UP \ x)
    ldiv!(sqrt(d.κ), x)
    return x
end

function Distributions._logpdf(d::GMRF, x::AbstractVector{T}) where T<:Real
    n = length(d)
    chol = cholesky(d.S)
    logpdf = - 0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * (n * log(d.κ) + logdet(chol))
    logpdf -= 0.5 * d.κ * x' * d.S * x
    return(logpdf)
end

# SUITE SPARSE FIXES TO WORK WITH SUBARRAYS -------------------------------------

# This functions is required for the default Distributions.rand!
function (\)(L::FactorComponent, b::AbstractVector)
    reshape(Matrix(L \ Dense(b)), length(b))
end

# PENALTY STRUCTURES BASES ------------------------------------------------------

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

# PENALTY STRUCTURES ------------------------------------------------------------

function structure(g::CartesianGrid{1}; δ = 0, order = 1) where {T}
    A = adjacency(g)
    if order == 1
        (2 + δ)I - A
    else
        throw(ErrorException("not implemented"))
    end
end

function structure_cyclic(g::CartesianGrid{1}; δ = 0, order = 1) where {T}
    A = adjacency_cyclic(g)
    if order == 1
        (2 + δ)I - A
    else
        throw(ErrorException("not implemented"))
    end
end
