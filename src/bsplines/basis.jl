"""
    basis(x, bs::RegularBsplines)

Returns an sparse evaluation of the basis splines `bs` at `x`.

The returning object is sparse because at any value, there are only `bs.order` non-zero
basis splines. The output is a sparse vector of size `bs.df` if `x` is a `Number`, and is
a sparse matrix of size `length(x)×bs.df` if `x` is a `Vector` or `AbstractRange`.
"""
function basis(x::Union{Number,AbstractRange,Vector}, bs::RegularBsplines)
    # 1) computes the `bs.order` non-zero basis splines and the indices for the last
    # non-zero basis splines of each element in `x`.
    basis, lastindices = nonzerobasis(x, bs)
    # 2) converts the basis to a sparse array associated to the `bs.df` basis splines at
    # `x`.
    sparsebasis(basis, lastindices, bs.df)
end

"""
    nonzerobasis(x::Number, bs::RegularBsplines)

Returns the `bs.order` non-zero basis splines at `x` and the index of the last non-zero
basis spline.

The returning vector has `bs.order` elements corresponding to the non-zero basis splines at
`x`. The index of the last non-zero basis spline is also returned to know which basis
splines are non-zero. The algorithm is based on De Boor (2001, page 110).
"""
function nonzerobasis(x::Number, bs::RegularBsplines)
    # knots to compute basis splines
    knots = extendedknots(bs)
    knotrange = range(knots)

    # initialize basis
    basis = zeros(bs.order)

    # check x inside bsplines domain
    if !(bs.lower ≤ x ≤ bs.upper)
        throw(DomainError(x, "The values of x should lie in the range of bs."))
    end

    # s index such as knots[s] ≤ x < knots[s+1]
    s_x = get_index(x, knotrange)

    # order = 1
    basis[1] = 1

    # order > 1
    δᵣ = zeros(bs.order - 1)
    δₗ = zeros(bs.order - 1)

    # uses basis of order k to compute basis of order k+1
    for k = 1:(bs.order - 1)
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
    nonzerobasis(x::Union{AbstractRange,Vector}, bs::RegularBsplines)

Returns the `bs.order` non-zero basis splines for each element of `x` and indices for the
last non-zero basis splines.

The returning matrix has `length(x)` rows and `bs.order` columns, each row correspond to
the non-zero basis splines of each element in `x`. The indices of the last non-zero basis
splines are also returned to know which basis splines are non-zero. The algorithm is based
on De Boor (2001, page 110).
"""
function nonzerobasis(x::Union{AbstractRange,Vector}, bs::RegularBsplines)
    n = length(x)

    # knots to compute basis splines
    knots = extendedknots(bs)
    knotrange = range(knots)

    # initialize basis (n × bs.order) and indices (s: knots[s] ≤ x < knots[s+1])
    basis = zeros(n, bs.order)
    indices = zeros(Int, n)

    # check x inside bsplines domain
    if !all(bs.lower ≤ val ≤ bs.upper for val in x)
        throw(DomainError(x, "The values of x should lie in the range of bs."))
    end

    # computing 'bs.order' non-zero functions at each 'x'
    for i in 1:n
        indices[i] = get_index(x[i], knotrange)
        s_x = indices[i]

        # order = 1
        basis[i, 1] = 1

        # order > 1
        δᵣ = zeros(bs.order - 1)::Vector
        δₗ = zeros(bs.order - 1)

        # uses basis of order k to compute basis of order k+1
        for k = 1:(bs.order - 1)
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
    sparsebasis(basis, lastindices, df::Int)

Converts the non-zero `basis` splines to an sparse array using the `lastindices` of the
non-zero basis.

The returning object is a sparse vector of dimension `df` if `basis` is a `Vector` and a
sparse matrix of dimension `size(basis, 1) × df` if the `basis` is a `Matrix`.

Given that the last non-zero basis spline has index `s`, it is known that the non-zero
function are those with indices (s, s-1, ..., s-order+1).
"""
function sparsebasis(basis::Vector, lastindex::Int, df::Int)
    order = length(basis)                      # number of non-zero basis splines
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

"""
    get_index(x::Number, knots::AbstractRange)

Get index `s` such as `knots[s] ≤ x < knots[s+1]`. Note that the basis spline s-th is the
last function that is non-zero. We use it to identify basis splines of certain `order`
that are non-zero at x.
"""
function get_index(x::Number, knots::AbstractRange)

    # check if x is inside the range
    if x < knots[1]
        return 0
    elseif x >= knots[end]
        return length(knots)
    end

    # compute index if x is inside the range
    s = 1 + floor(Int, (x - knots[1]) / step(knots))

    # simple fix in case condition knots[s] ≤ x < knots[s+1] is not hold
    if x >= knots[s + 1]
        s = s + 1
    end
    if x < knots[s]
        s = s - 1
    end

    return s
end

# ------------------------------
# Old functions. Needs cleaning.
# ------------------------------

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
    j = get_index(x, knotrange)
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
        j_x = get_index(x[i], knotrange)
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
