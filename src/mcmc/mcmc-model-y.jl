function sample_model_y(y, Bw, Pw, κw, id; binit = nothing, σ²yinit = nothing, δwinit = nothing,
    nsamples = 10, burnin = 0, thin = 1, thinz = 1)

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

    niter = mcmcniter(nsamples, burnin, thin)
    nsamplesz = mcmcnsamples(niter, burnin, thinz)

    # initialize samples
    δw_samples = zeros(qw, nsamples + 1)
    z_samples = [zeros(n[k], nsamplesz + 1) for k in 1:K]
    b_samples = zeros(K, nsamples + 1)
    σ²y_samples = zeros(K, nsamples + 1)

    δw_samples[:, 1] = δw
    [z_samples[k][:, 1] = z[k] for k in 1:K]
    b_samples[:, 1] = b
    σ²y_samples[:, 1] = σ²y

    # mcmc
    for i = 1:niter

        # check if iteration should be saved
        saveiter = mcmcsave(i, burnin, thin)
        saveid = mcmcid(i, burnin, thin) + 1
        saveiterz = mcmcsave(i, burnin, thinz)
        saveidz = mcmcid(i, burnin, thinz) + 1

        # sample z
        [resid[k] -= z[k] for k in 1:K]
        μz = [-resid[k] for k in 1:K]
        varz = σ²y
        Z = [truncated.(Normal.(μz[k], sqrt.(varz[k])), zlower[k], zupper[k]) for k in 1:K]
        z  = [rand.(Z[k]) for k in 1:K]
        if saveiterz
            println("Iteration $i with saved sample $(saveid):")
            for k = 1:K
                z_samples[k][:, saveidz] = z[k]
            end
        end
        [resid[k] += z[k] for k in 1:K]

        # sample δw
        Qw = sum(BwtBw ./ σ²y) + κw * Pw
        CQw = cholesky(Qw)
        [resid[k] += Bw[k] * δw for k in 1:K]
        aux = sum([σ²y[k]^(-1) * Bw[k]' * resid[k] for k in 1:K])
        μw =  CQw.UP \ (CQw.PtL \ aux)
        δw = μw + CQw.UP \ randn(qw)
        if saveiter
            δw_samples[:, saveid] = δw
        end
        [resid[k] -= Bw[k] * δw for k in 1:K]

        # sample β
        Qᵦ = sum([σ²y[k]^(-1) * A[k]' * A[k] for k in 1:K]) + Σᵦ
        CQᵦ = cholesky(Qᵦ)
        μᵦ = CQᵦ.U \ (CQᵦ.L \ sum([σ²y[k]^(-1) * A[k]' * (z[k] - Bw[k] * δw) for k in 1:K]))
        b = μᵦ + CQᵦ.U \ randn(K)
        if saveiter
            b_samples[:, saveid] = b
        end
        resid = [z[k] - A[k] * b - Bw[k] * δw for k in 1:K]

        # sample σ²y
        a_σ²y = sigma2_a_prior .+ n / 2
        b_σ²y = [sigma2_b_prior + sum(resid[k].^2) / 2 for k in 1:K]
        σ²y = rand.(InverseGamma.(a_σ²y, b_σ²y))
        if saveiter
            σ²y_samples[:, saveid] = σ²y
        end
    end

    Dict("w" => δw_samples, "beta" => b_samples, "sigma2_y" => σ²y_samples, "z" => z_samples)
end

function sample_model_y_fixsigma(y, Bw, Pw, κw, id; binit = nothing, σ²yinit = nothing, δwinit = nothing,
    κwinit = nothing, nsamples = 10, burnin = 0, thin = 1, thinz = 1)

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
    κw = isnothing(κwinit) ? 1.0 : copy(κwinit)

    # pre-computation
    BwtBw = [Bw[k]' * Bw[k] for k in 1:K]

    # using design matrix approach
    A = [myA(k, n[k], K, id) for k in 1:K]
    resid = [z[k] - A[k] * b - Bw[k] * δw for k in 1:K]

    # truncated limits
    zlower = [ifelse.(y[k] .== 0, nothing, 0) for k in 1:K]
    zupper = [ifelse.(y[k] .== 1, nothing, 0) for k in 1:K]

    niter = mcmcniter(nsamples, burnin, thin)
    nsamplesz = mcmcnsamples(niter, burnin, thinz)

    # initialize samples
    δw_samples = zeros(qw, nsamples + 1)
    z_samples = [zeros(n[k], nsamplesz + 1) for k in 1:K]
    b_samples = zeros(K, nsamples + 1)
    σ²y_samples = zeros(K, nsamples + 1)
    κw_samples = zeros(1, nsamples + 1)

    δw_samples[:, 1] = δw
    [z_samples[k][:, 1] = z[k] for k in 1:K]
    b_samples[:, 1] = b
    σ²y_samples[:, 1] = σ²y
    κw_samples[1, 1] = κw

    # mcmc
    for i = 1:niter

        # check if iteration should be saved
        saveiter = mcmcsave(i, burnin, thin)
        saveid = mcmcid(i, burnin, thin) + 1
        saveiterz = mcmcsave(i, burnin, thinz)
        saveidz = mcmcid(i, burnin, thinz) + 1

        # sample z
        [resid[k] -= z[k] for k in 1:K]
        μz = [-resid[k] for k in 1:K]
        varz = σ²y
        Z = [truncated.(Normal.(μz[k], sqrt.(varz[k])), zlower[k], zupper[k]) for k in 1:K]
        z  = [rand.(Z[k]) for k in 1:K]
        if saveiterz
            println("Iteration $i with saved sample $(saveid):")
            for k = 1:K
                z_samples[k][:, saveidz] = z[k]
            end
        end
        [resid[k] += z[k] for k in 1:K]

        # sample δw
        Qw = sum(BwtBw ./ σ²y) + κw * Pw
        CQw = cholesky(Qw)
        [resid[k] += Bw[k] * δw for k in 1:K]
        aux = sum([σ²y[k]^(-1) * Bw[k]' * resid[k] for k in 1:K])
        μw =  CQw.UP \ (CQw.PtL \ aux)
        δw = μw + CQw.UP \ randn(qw)
        if saveiter
            δw_samples[:, saveid] = δw
        end
        [resid[k] -= Bw[k] * δw for k in 1:K]

        # sample β
        Qᵦ = sum([σ²y[k]^(-1) * A[k]' * A[k] for k in 1:K]) + Σᵦ
        CQᵦ = cholesky(Qᵦ)
        μᵦ = CQᵦ.U \ (CQᵦ.L \ sum([σ²y[k]^(-1) * A[k]' * (z[k] - Bw[k] * δw) for k in 1:K]))
        b = μᵦ + CQᵦ.U \ randn(K)
        if saveiter
            b_samples[:, saveid] = b
        end
        resid = [z[k] - A[k] * b - Bw[k] * δw for k in 1:K]

        # sample σ²y
        if K > 1
            a_σ²y = sigma2_a_prior .+ n[2:K] / 2
            b_σ²y = [sigma2_b_prior + sum(resid[k].^2) / 2 for k in 2:K]
            σ²y[2:K] = rand.(InverseGamma.(a_σ²y, b_σ²y))
        end
        if saveiter
            σ²y_samples[:, saveid] = σ²y
        end

        # sample κw
        a_κ = kappa_a_prior + qw / 2
        b_κ = kappa_b_prior + (transpose(δw) * Pw * δw) / 2
        κw = rand(Gamma(a_κ, 1 / b_κ))
        if saveiter
            κw_samples[1, saveid] = κw
        end
    end

    Dict("w" => δw_samples, "beta" => b_samples, "sigma2_y" => σ²y_samples, "z" => z_samples, "κw" => κw_samples)
end

