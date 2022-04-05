"""
    RGMRF(grid, order, δ, κ)

Regular Gaussian Markov random field with n-th order, δ as and κ as precision.
"""
struct RGMRF <: AbstractGMRF
    grid::CartesianGrid
    order::Integer
    δ::Number
    κ::Real
end

Base.length(d::RGMRF) = nelements(d.grid)
structure(d::RGMRF) = structure(d.grid; δ = d.δ, order = d.order, cyclic = false)
scale(d::RGMRF) = d.κ
