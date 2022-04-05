"""
    GMRF(S, κ)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: AbstractGMRF
    S::AbstractMatrix
    κ::Real
end

Base.length(d::GMRF) = size(d.S, 1)
structure(d::GMRF) = d.S
scale(d::GMRF) = d.κ
