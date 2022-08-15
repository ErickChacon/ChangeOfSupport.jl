"""
    basis(x, bs::RegularBsplines)

Return an sparse vector of the evaluation of the basis splines `bs` at `x`.

The output is a sparse vector of size `bs.df` is x is a `Number`, and is a sparse matrix
of size `n×bs.df` if x is a `Vector`. The algorithm is based on De Boor (2001, page 110).
"""
function basis(x::Number, bs::RegularBsplines)
    # 1) we compute only the `bs.order` non-zero basis splines and the index for the last
    # nons-zero basis spline.
    basis, lastindices = nonzerobasis(x, bs)
    # 2) using the non-zero basis splines and the index, we create an sparse vector
    # associated to the `bs.df` basis splines.
    sparsebasis(basis, lastindices, bs.df)
end

"""
    basis(x::AbstractVector, bs::RegularBsplines)

Return an sparse matrix of the evaluation of Bsplines `bs` at `x`.

The sparse matrix has `length(x)` rows and `bs.df` columns. The algorithm is based on De
Boor (2001, page 110). There are only `bs.order` non-zero functions at each value of `x`.
"""
function basis(x::Union{AbstractRange,Vector}, bs::RegularBsplines)
    # 1) we compute only the `bs.order` non-zero basis splines and the indices for the
    # last nons-zero basis spline of each element in `x`.
    basis, lastindices = nonzerobasis(x, bs)
    # 2) we convert the basis to a sparse design matrix associated to the `bs.df` basis
    # splines.
    sparsebasis(basis, lastindices, bs.df)
end

"""
    nonzerobasis(x::Number, b::RegularBsplines)

Returns the `b.order` non-zero basis splines at `x` and the index of the last non-zero
basis splines.

The returning vector has `b.order` elements corresponding to the non-zero basis splines at
`x`. The algorithm is based on De Boor (2001, page 110).
"""
function nonzerobasis(x::Number, b::RegularBsplines)
    # knots to compute basis splines
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis
    basis = zeros(b.order)

    # check x inside bsplines domain
    if !(b.lower ≤ x ≤ b.upper)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # s index such as knots[s] ≤ x < knots[s+1]
    s_x = get_x_index(x, knotrange)

    # order = 1
    basis[1] = 1

    # order > 1
    δᵣ = zeros(b.order - 1)
    δₗ = zeros(b.order - 1)

    # uses basis of order k to compute basis of order k+1
    for k = 1:(b.order - 1)
        # create left and right differences
        δᵣ[k] = knotrange[s_x + k] - x
        δₗ[k] = x - knotrange[s_x + 1 - k]
        leftterm = 0
        # evaluate the first k non-zero basis of order k+1 (r = 1:k)
        for r = 1:k
            zs = basis[r] / (δᵣ[r] + δₗ[k-r+1])
            basis[r] = leftterm + δᵣ[r] * zs
            # leftterm for r' = r+1, then δ(k-r'+2) = δ(k-r+1) and z(r'-1) = z(r).
            leftterm = δₗ[k-r+1] * zs
        end
        # evaluate the k+1 basis of order k+1 (r = k+1)
        basis[k+1] = leftterm
    end

    # return the non-zero basis splines and index for last non-zero basis
    basis, s_x
end

"""
    nonzerobasis(x::Union{AbstractRange,Vector}, b::RegularBsplines)

Returns the `b.order` non-zero basis splines for each element of `x` and indices for the last non-zero basis splines.

The returning matrix has `length(x)` rows and `b.order` columns, each row correspond to
the non-zero basis splines of each element in `x`. The algorithm is based on De Boor
(2001, page 110).
"""
function nonzerobasis(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    n = length(x)

    # knots to compute basis splines
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis (n × b.order) and indices (s: knots[s] ≤ x < knots[s+1])
    basis = zeros(n, b.order)
    indices = zeros(Int, n)

    # check x inside bsplines domain
    if !all(b.lower ≤ val ≤ b.upper for val in x)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # computing 'b.order' non-zero functions at each 'x'
    for i in 1:n
        indices[i] = get_x_index(x[i], knotrange)
        s_x = indices[i]

        # order = 1
        basis[i, 1] = 1

        # order > 1
        δᵣ = zeros(b.order - 1)
        δₗ = zeros(b.order - 1)

        # uses basis of order k to compute basis of order k+1
        for k = 1:(b.order - 1)
            # create left and right differences
            δᵣ[k] = knotrange[s_x + k] - x[i]
            δₗ[k] = x[i] - knotrange[s_x + 1 - k]
            leftterm = 0
            # evaluate the first k non-zero basis of order k+1 (r = 1:k)
            for r = 1:k
                zs = basis[i, r] / (δᵣ[r] + δₗ[k-r+1])
                basis[i, r] = leftterm + δᵣ[r] * zs
                # leftterm for r' = r+1, then δ(k-r'+2) = δ(k-r+1) and z(r'-1) = z(r).
                leftterm = δₗ[k+1-r] * zs
            end
            # evaluate the kExact representaions willl demand a different abstract type that is not Integer.+1 basis of order k+1 (r = k+1)
            basis[i, k + 1] = leftterm
        end
    end

    # return the non-zero basis splines and index for last non-zero basis
    basis, indices
end

"""
    sparsebasis(basis::Vector, lastindex::Int, df::Int)

Converts the non-zero basis splines to an sparse matrix of dimension `n×df`.
"""
function sparsebasis(basis::Vector, lastindex::Int, df::Int)
    order = length(basis) # number of non-zero basis splines
    J = [(lastindex-order+k) for k in 1:order] # indices of basis splines
    sparsevec(J, basis, df + 1)[1:df]
end

function sparsebasis(basis::Matrix, lastindices::Vector, df::Int)
    (n, order) = size(basis)
    I = repeat(1:n, outer = order)                          # index of observations
    J = [(s-order+k) for k in 1:order for s in lastindices] # index of basis splines
    sparse(I, J, vec(basis), n, df + 1)[:, 1:df]
end

# we can improve sparsebasis considering the restriction (lower ≤ x < upper) instead of
# (lower ≤ x ≤ upper).

# The following functions are just old functions used to be compared with current
# versions.

function basis_sparse(x::Number, b::RegularBsplines)
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis
    basis = spzeros(length(knotrange))

    # check x inside acceptable range
    if !(b.lower ≤ x ≤ b.upper)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # order = 1
    j = get_x_index(x, knotrange)
    basis[j] = 1

    # order > 1
    δᵣ = zeros(b.order - 1)
    δₗ = zeros(b.order - 1)

    # uses basis of order k to compute basis of order k+1
    for k = 1:(b.order - 1)
        # create left and right differences
        δᵣ[k] = knotrange[j + k] - x
        δₗ[k] = x - knotrange[j + 1 - k]
        leftterm = 0
        # evaluate the first k non-zero basis of order k+1 (s = 1:k)
        for s = 1:k
            zs = basis[j-k+s] / (δᵣ[s] + δₗ[k-s+1])
            basis[j-k-1+s] = leftterm + δᵣ[s] * zs
            # leftterm for s' = s+1, then δ(k-s'+2) = δ(k-s+1) and z(s'-1) = z(s).
            leftterm = δₗ[k-s+1] * zs
        end
        # evaluate the k+1 basis of order k+1 (s = k+1)
        basis[j] = leftterm
    end
    # return only the first df basis functions, these are non-zero in the desired range
    basis[1:b.df]
end

function basis_sparse(x::Union{AbstractRange,Vector}, b::RegularBsplines)
    n = length(x)

    # knots to compute basis splines
    knots = extendedknots(b)
    knotrange = range(knots)

    # initialize basis and indices to create sparse arrays (n × b.df size)
    basis = zeros(n * b.order) # basis values
    I = repeat(1:n, inner = b.order) # index of observation
    J = zeros(Int, n * b.order) # index of basis function

    # check x inside bsplines domain
    if !all(b.lower ≤ val ≤ b.upper for val in x)
        throw(DomainError(x, "The values of x should lie in the range of b."))
    end

    # computing 'b.order' non-zero functions at each 'x'
    for i in 1:n
        # order = 1
        ibase = (i-1) * b.order
        j_x = get_x_index(x[i], knotrange)
        basis[ibase + 1] = 1

        # order > 1
        δᵣ = zeros(b.order - 1)
        δₗ = zeros(b.order - 1)

        # uses basis of order k to compute basis of order k+1
        for k = 1:(b.order - 1)
            # create left and right differences
            δᵣ[k] = knotrange[j_x + k] - x[i]
            δₗ[k] = x[i] - knotrange[j_x + 1 - k]
            leftterm = 0
            # evaluate the first k non-zero basis of order k+1 (s = 1:k)
            for s = 1:k
                zs = basis[ibase + s] / (δᵣ[s] + δₗ[k+1-s])
                basis[ibase + s] = leftterm + δᵣ[s] * zs
                # leftterm for s' = s+1, then δ(k-s'+2) = δ(k-s+1) and z(s'-1) = z(s).
                leftterm = δₗ[k+1-s] * zs
            end
            # evaluate the k+1 basis of order k+1 (s = k+1)
            basis[ibase + k + 1] = leftterm
        end

        # J indices for non-zero basis
        J[(ibase+1):(ibase+b.order)] = (j_x-b.order+1):j_x
    end

    # return only the first df basis functions
    sparse(I, J, basis, length(x), length(knotrange))[:, 1:b.df]
end
