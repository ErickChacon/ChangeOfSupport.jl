
function sample_gam(y, X, P, σ², κ, niter = 10; fix_hyper = false)
    n, p = size(X)
    XtX = X' * X

    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001

    β_samples = zeros(p, niter)
    σ²_samples = zeros(1, niter)
    κ_samples = zeros(1, niter)

    σ²_samples[1, 1] = σ²
    κ_samples[1, 1] = κ

    # out = 1

    for i = 2:niter
        # sample β
        Qᵦ = XtX / σ² + κ * P
        CQᵦ = cholesky(Qᵦ)
        μᵦ =  CQᵦ.U \ (CQᵦ.L \ (transpose(X) * y / σ²))
        β = μᵦ + CQᵦ.U \ randn(p)
        β_samples[:, i] = β

        # sample σ²
        a_σ² = sigma2_a_prior + n / 2
        b_σ² = sigma2_b_prior + sum((y - X * β).^2) / 2
        σ² = rand(InverseGamma(a_σ², b_σ²))
        σ²_samples[1, i] = σ²

        if !fix_hyper
            # sample κ
            a_κ = kappa_a_prior + (p - 1) / 2
            b_κ = kappa_b_prior + transpose(β) * P * β / 2
            κ = rand(Gamma(a_κ, 1 / b_κ))
        end
        κ_samples[1, i] = κ
    end

    β_samples, σ²_samples, κ_samples
    # out
end

function sample_gam_sparse(y, X, P, σ², κ, niter = 10; fix_hyper = false)
    n, p = size(X)
    XtX = X' * X

    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001

    β_samples = zeros(p, niter)
    σ²_samples = zeros(1, niter)
    κ_samples = zeros(1, niter)

    σ²_samples[1, 1] = σ²
    κ_samples[1, 1] = κ

    # out = 1

    for i = 2:niter
        # sample β
        Qᵦ = XtX / σ² + κ * P
        CQᵦ = cholesky(Qᵦ)
        μᵦ =  CQᵦ.UP \ (CQᵦ.PtL \ (X'y / σ²))
        β = μᵦ + CQᵦ.UP \ randn(p)
        β_samples[:, i] = β

        # sample σ²
        a_σ² = sigma2_a_prior + n / 2
        b_σ² = sigma2_b_prior + sum((y - X * β).^2) / 2
        σ² = rand(InverseGamma(a_σ², b_σ²))
        σ²_samples[1, i] = σ²

        if !fix_hyper
            # sample κ
            a_κ = kappa_a_prior + (p - 1) / 2
            b_κ = kappa_b_prior + transpose(β) * P * β / 2
            κ = rand(Gamma(a_κ, 1 / b_κ))
        end
        κ_samples[1, i] = κ
    end

    β_samples, σ²_samples, κ_samples
    # out
end

function sample_gam_sparse2(y::AbstractVector{Float64}, X::SparseMatrixCSC{Float64, Int64},
    P::SparseMatrixCSC{Int64, Int64}, σ²::Float64, κ::Float64, niter::Int64 = 10; fix_hyper = false)
    n, p = size(X)
    XtX = X' * X

    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001

    β_samples = zeros(p, niter)
    σ²_samples = zeros(1, niter)
    κ_samples = zeros(1, niter)

    σ²_samples[1, 1] = σ²
    κ_samples[1, 1] = κ

    auxchol = cholesky(XtX / σ² + κ * P)
    perm_index = auxchol.p

    # out = 1

    for i = 2:niter
        # sample β
        @timeit "Qb" Qᵦ = XtX / σ² + κ * P
        @timeit "cholesky" CQᵦ = cholesky(Qᵦ, perm = perm_index)
        @timeit "mu" μᵦ =  CQᵦ.UP \ (CQᵦ.PtL \ (X'y / σ²))
        @timeit "beta" β = μᵦ + CQᵦ.UP \ randn(p)
        @timeit "store beta" β_samples[:, i] = β

        # sample σ²
        @timeit "a" a_σ² = sigma2_a_prior + n / 2
        @timeit "b" b_σ² = sigma2_b_prior + sum((y - X * β).^2) / 2
        @timeit "sigma^2" σ² = rand(InverseGamma(a_σ², b_σ²))
        @timeit "store sigma^2" σ²_samples[1, i] = σ²

        if !fix_hyper
            # sample κ
            @timeit "a_k" a_κ = kappa_a_prior + (p - 1) / 2
            @timeit "b_k" b_κ = kappa_b_prior + transpose(β) * P * β / 2
            @timeit "kappa" κ = rand(Gamma(a_κ, 1 / b_κ))
        end
        @timeit "store kappa" κ_samples[1, i] = κ
    end

    print_timer()
    β_samples, σ²_samples, κ_samples
    # out
end


function sample_gam_area(y, X, P, area, σ², κ, niter = 10; fix_hyper = false)
    n, p = size(X)
    Di_sqrt = Diagonal(sqrt.(area))
    Di_sqrtX = Di_sqrt * X
    XDiX = transpose(Di_sqrtX) * Di_sqrtX

    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001

    β_samples = zeros(p, niter)
    σ²_samples = zeros(1, niter)
    κ_samples = zeros(1, niter)

    σ²_samples[1, 1] = σ²
    κ_samples[1, 1] = κ

    for i = 2:niter
        # sample β
        Qᵦ = XDiX / σ² + κ * P
        CQᵦ = cholesky(Qᵦ)
        μᵦ =  CQᵦ.U \ (CQᵦ.L \ (transpose(X) * (area .* y) / σ²))
        β = μᵦ + CQᵦ.U \ randn(p)
        β_samples[:, i] = β

        # sample σ²
        a_σ² = sigma2_a_prior + n / 2
        b_σ² = sigma2_b_prior + sum(area .* (y - X * β) .^ 2) / 2
        σ² = rand(InverseGamma(a_σ², b_σ²))
        σ²_samples[1, i] = σ²

        if !fix_hyper
            # sample κ
            a_κ = kappa_a_prior + (p - 1) / 2
            b_κ = kappa_b_prior + transpose(β) * P * β / 2
            κ = rand(Gamma(a_κ, 1 / b_κ))
        end
        κ_samples[1, i] = κ
    end

    β_samples, σ²_samples, κ_samples
end
