"""
    CGMRF(grid, order, δ, κ)

Circulant Gaussian Markov random fields with n-th order, δ as and κ as precision.
"""
struct CGMRF <: AbstractGMRF
    grid::CartesianGrid
    order::Integer
    δ::Number
    κ::Number
end

Base.length(d::CGMRF) = nelements(d.grid)
structure(d::CGMRF) = structure(d.grid; δ = d.δ, order = d.order, cyclic = true)
scale(d::CGMRF) = d.κ
structure_base(d::CGMRF) = structure_base(d.grid; δ = d.δ, order = d.order)


## Random generator

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractVector{T}) where T<:Real
    # get dimensions
    dims = size(d.grid)
    n = prod(dims)

    # compute eigenvalues
    base = structure_base(d)
    λ = FFTW.fft(base)

    # simulate using fft
    z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
    xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
    if length(dims) == 2
        xaux = transpose(reverse(xaux, dims = 1))
    end
    copyto!(x, xaux)
    ldiv!(sqrt(n * scale(d)), x)

    return x
end

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractArray{<:Real})
    dims = size(d.grid)
    n = prod(dims)
    λ = FFTW.fft(structure_base(d))
    @inbounds for xi in Distributions.eachvariate(x, Distributions.variate_form(typeof(d)))
        z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
        xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
        if length(dims) == 2
            xaux = transpose(reverse(xaux, dims = 1))
        end
        copyto!(xi, xaux)
        ldiv!(sqrt(n * scale(d)), xi)
    end
    return x
end

# Logpdf

# function Distributions._logpdf(X::CGMRF, x::AbstractVector{T}) where T<:Real
#     # get dimensions
#     dims = size(X.grid)
#     n = prod(dims)
#
#     # compute eigenvalues
#     base = structure_base(X.grid; δ = X.δ, order = X.order)
#     λ = fft(base)
#
#     # output = (- n * log(2π) + sum(log.(real(λ)))) / 2
#     output = λ .* ifft()
# end
