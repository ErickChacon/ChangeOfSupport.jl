# This file needs to be modified to represent a general IGMRF.

"""
    IGMRF(grid, order)

Regular Gaussian Markov random fields.
"""
struct IGMRF{N,T}
    grid::CartesianGrid{N,T}
    order::Number
end

# function to convert i,j cel to k-index
function ij_to_k(i, j, n1, n2)
    j + (i - 1) * n2
end

# function to obtain the precision matrix of a igmrf
function igmrf_precision_1(n1; δ = 0)

    n = n1

    # off-diagonal
    I1 = 1:n1-1
    J1 = 2:n
    I2 = 2:n1
    J2 = 1:n1-1
    D1 = [1, n1]
    D2 = 2:n1-1

    # sparse precision matrix up to a constant
    P = sparse(
        vcat(I1, I2, D1, D2),
        vcat(J1, J2, D1, D2),
        vcat(-1 * ones(2n1-2), ones(2), (2 + δ) * ones(n1-2)),
        n1, n1
     )

    return P
end

# function to obtain the precision matrix of a igmrf
function igmrf_precision_1_t(n1; δ = 0)

    n = n1

    # off-diagonal
    I1 = 1:n1-1
    J1 = 2:n
    I2 = 2:n1
    J2 = 1:n1-1
    D1 = [1, n1]
    D2 = 2:n1-1

    # sparse precision matrix up to a constant
    P = sparse(
        vcat(I1, I2, D1, D2),
        vcat(J1, J2, D1, D2),
        vcat(-1 * ones(2n1-2), (2 + δ) * ones(2), (2 + δ) * ones(n1-2)),
        n1, n1
     )

    return P
end


# function to obtain the precision matrix of a igmrf
function igmrf_precision_1(n1, n2)

    n = n1 * n2

    # off-diagonal
    I1 = [ij_to_k(i, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
    J1 = [ij_to_k(i + 1, j, n1, n2) for i = 1:(n1-1) for j = 1:n2]
    I2 = [ij_to_k(i, j, n1, n2) for i = 2:n1 for j = 1:n2]
    J2 = [ij_to_k(i - 1, j, n1, n2) for i = 2:n1 for j = 1:n2]
    I3 = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
    J3 = [ij_to_k(i, j + 1, n1, n2) for i = 1:n1 for j = 1:(n2-1)]
    I4 = [ij_to_k(i, j, n1, n2) for i = 1:n1 for j = 2:n2]
    J4 = [ij_to_k(i, j - 1, n1, n2) for i = 1:n1 for j = 2:n2]

    # diagonal
    D2 = [ij_to_k(i, j, n1, n2) for (i,j) in [(1,1), (1,n2), (n1,1), (n1,n2)]]
    D3_1 = [ij_to_k(i, j, n1, n2) for i in [1, n1] for j in 2:(n2-1)]
    D3_2 = [ij_to_k(i, j, n1, n2) for j in [1, n2] for i in 2:(n1-1)]
    D4 = [ij_to_k(i, j, n1, n2) for i = 2:(n1-1) for j = 2:(n2-1)]

    # sparse precision matrix up to a constant
    n_off = 2(2 * n1 * n2 - n1 - n2)
    n_d2 = 4
    n_d3 = 2(n1 + n2 - 4)
    n_d4 = n1 * n2 - 2n1 - 2n2 + 4
    P = sparse(
        vcat(I1, I2, I3, I4, D2, D3_1, D3_2, D4),
        vcat(J1, J2, J3, J4, D2, D3_1, D3_2, D4),
        vcat(-1 * ones(n_off), 2 * ones(n_d2), 3 * ones(n_d3), 4 * ones(n_d4)),
        n, n
     )

    return P
end

# build circulant base of a gmrf
function igmrf_precision_base(n::Int; δ = 0.01, k = 1)
    base = zeros(n)
    if k == 1
        base[[1, 2, end]] = [2 + δ, -1, -1]
    elseif k == 2
        base[[1, 2, 3, end-1, end]] = [6 + δ, -4, 1, 1, -4]
    end
    return base
end

# build circulant base of a gmrf
function igmrf_precision_base(n1::Int, n2::Int; δ = 0.01, k = 1)
    n = n1 * n2
    base = zeros(n2, n1)
    if k == 1
        base[1, 1] = 4 + δ
        base[1, 2] = -1
        base[2, 1] = -1
        base[1, end] = -1
        base[end, 1] = -1
    elseif k == 2
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
    end
    return base
end

# simulate igmrf
function igmrf_simulate(dims; δ = 0.01, k = 1, border = 0)
    n = prod(dims)

    # define base
    base = igmrf_precision_base(dims...; δ = δ, k = k)

    # compute eigenvalues
    λ = fft(base)

    # simulate
    z = randn(reverse(dims)...) + randn(reverse(dims)...) * im
    x = real(fft(λ .^ (-1/2) .* z) / sqrt(n))
    #
    # # edge effects
    # if length(dims) == 1
    #     return x[border + 1:n - border]
    # elseif length(dims) == 2
    #     return x[border + 1:dims[2] - border, border + 1:dims[1] - border]
    # end
end

function simulate(igmrf::IGMRF{N,T}; δ = 0.01, border = 0) where {N,T}
    dims = igmrf.grid.dims
    n = prod(dims)

    # define base
    base = igmrf_precision_base(dims..., k = igmrf.order)

    # compute eigenvalues
    λ = sqrt(n) * fft(base)

    # simulate
    z = randn(reverse(dims)...) + randn(reverse(dims)...) * im
    x = real(fft(λ .^ (-1/2) .* z))

    # edge effects
    if length(dims) == 1
        return x[border + 1:n - border]
    elseif length(dims) == 2
        return x[border + 1:dims[2] - border, border + 1:dims[1] - border]
    end
end
