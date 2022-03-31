"""
    CGMRF(grid, order, δ, κ)

Circulant Gaussian Markov random fields with n-th order, δ as and κ as precision.
"""
# Sampleable{Multivariate, Continuous}
struct CGMRF{N,T} <: ContinuousMultivariateDistribution
    grid::CartesianGrid{N,T}
    order::Integer
    δ::Number
    κ::Number
end


function structure(X::CGMRF{1,T}) where {T}
    A = adjacency(X.grid)
    if X.order == 1
        (2 + X.δ)I - A
    end
end

Base.length(X::CGMRF) = nelements(X.grid)

# function Distributions._rand!(rng::AbstractRNG, X::CGMRF, x::AbstractVector{T}) where T<:Real
#     # get dimensions
#     dims = size(X.grid)
#     n = prod(dims)
#
#     # compute eigenvalues
#     base = structure_base(X.grid; δ = X.δ, order = X.order)
#     λ = fft(base)
#
#     # simulate using fft
#     z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
#     copyto!(x, real(fft(λ .^ (-1/2) .* z)))
#     ldiv!(sqrt(n * X.κ), x)
#
#     return x
# end

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

# function Distributions._logpdf(X::CGMRF, x::AbstractVector{T}) where T<:Real
#
# end

# function Distributions._logpdf(X::CGMRF, x::AbstractVector{T}) where T<:Real
#     # # get dimensions
#     # dims = size(X.grid)
#     # n = prod(dims)
#     #
#     # # compute eigenvalues
#     # base = structure_base(X.grid; δ = X.δ, order = X.order)
#     # λ = fft(base)
#     #
#     # # output = (- n * log(2π) + sum(log.(real(λ)))) / 2
#     # output = λ .* ifft()
#     #
#     # # shdfhdim, v = mvtdist_consts(d)
#     # # v - shdfhdim * log1p(sqmahal(d, x) / d.df)
# end

# function _logpdf!(r::AbstractArray{<:Real}, d::AbstractMvNormal, x::AbstractMatrix{<:Real})
#     sqmahal!(r, d, x)
#     c0 = mvnormal_c0(d)
#     for i = 1:size(x, 2)
#         @inbounds r[i] = c0 - r[i]/2
#     end
#     r
# end
