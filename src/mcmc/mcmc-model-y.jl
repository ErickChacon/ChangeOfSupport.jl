function sample_model_y(y, Bw, Pw, κw, id; binit = nothing, σ²yinit = nothing, δwinit = nothing, niter = 10)

    # dimensions
    K = length(y)
    n = length.(y)
    qw = size(Bw[1], 2)

    # hyperparameters
    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001
    Σᵦ = 3I

    # initial values
    z = [zeros(n[k]) for k in 1:K]
    b = isnothing(binit) ? zeros(K) : copy(binit)
    σ²y = isnothing(σ²yinit) ? ones(K) : copy(σ²yinit)
    δw = isnothing(δwinit) ? zeros(qw) : copy(δwinit)

    # pre-computation
    BwtBw = [Bw[k]' * Bw[k] for k in 1:K]

    # using design matrix approach
    A = [myA(k, n[k], K, id) for k in 1:K]
    resid = [z[k] - A[k] * b - Bw[k] * δw for k in 1:K]

    # truncated limits
    zlower = [ifelse.(y[k] .== 0, nothing, 0) for k in 1:K]
    zupper = [ifelse.(y[k] .== 1, nothing, 0) for k in 1:K]

    # initialize samples
    δw_samples = zeros(qw, niter)
    z_samples = [zeros(n[k], niter) for k in 1:K]
    b_samples = zeros(K, niter)
    σ²y_samples = zeros(K, niter)

    δw_samples[:, 1] = δw
    [z_samples[k][:, 1] = z[k] for k in 1:K]
    b_samples[:, 1] = b
    σ²y_samples[:, 1] = σ²y

    # mcmc
    for i = 2:niter
        # sample z
        [resid[k] -= z[k] for k in 1:K]
        μz = [-resid[k] for k in 1:K]
        varz = σ²y
        Z = [truncated.(Normal.(μz[k], sqrt.(varz[k])), zlower[k], zupper[k]) for k in 1:K]
        z  = [rand.(Z[k]) for k in 1:K]
        for k = 1:K
            z_samples[k][:, i] = z[k]
        end
        [resid[k] += z[k] for k in 1:K]

        # sample δw
        Qw = sum(BwtBw ./ σ²y) + κw * Pw
        CQw = cholesky(Qw)
        [resid[k] += Bw[k] * δw for k in 1:K]
        aux = sum([σ²y[k]^(-1) * Bw[k]' * resid[k] for k in 1:K])
        μw =  CQw.UP \ (CQw.PtL \ aux)
        δw = μw + CQw.UP \ randn(qw)
        δw_samples[:, i] = δw
        [resid[k] -= Bw[k] * δw for k in 1:K]

        # sample β
        Qᵦ = sum([σ²y[k]^(-1) * A[k]' * A[k] for k in 1:K]) + Σᵦ
        CQᵦ = cholesky(Qᵦ)
        μᵦ = CQᵦ.U \ (CQᵦ.L \ sum([σ²y[k]^(-1) * A[k]' * (z[k] - Bw[k] * δw) for k in 1:K]))
        b = μᵦ + CQᵦ.U \ randn(K)
        b_samples[:, i] = b
        resid = [z[k] - A[k] * b - Bw[k] * δw for k in 1:K]

        # sample σ²y
        a_σ²y = sigma2_a_prior .+ n / 2
        b_σ²y = [sigma2_b_prior + sum(resid[k].^2) / 2 for k in 1:K]
        σ²y = rand.(InverseGamma.(a_σ²y, b_σ²y))
        σ²y_samples[:, i] = σ²y
    end

    δw_samples, b_samples, σ²y_samples, z_samples
end

