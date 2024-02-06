function sample_model_standard(y, x, Bw, Bvx, Bvy, Pw, Pv, κw, id; binit = nothing, βinit = nothing,
    σ²yinit = nothing, σ²xinit = nothing, κvinit = nothing, δwinit = nothing, δvinit = nothing,
    nsamples = 10, burnin = 0, thin = 1, thinz = 1)

    # dimensions
    K = length(y)
    p = length(x)
    pf = K + p
    n = length.(y)
    m = map(x -> size(x, 1), x)
    qw = size(Bw[1], 2)
    qv = map(x -> size(x, 2), Bvx)

    # hyperparameters
    sigma2_a_prior = 0.001
    sigma2_b_prior = 0.001
    kappa_a_prior = 0.001
    kappa_b_prior = 0.001
    Σᵦ = 3I

    # initial values
    z = [zeros(n[k]) for k in 1:K]
    b = isnothing(binit) ? zeros(K) : copy(binit)
    β = isnothing(βinit) ? zeros(p) : copy(βinit)
    σ²y = isnothing(σ²yinit) ? ones(K) : copy(σ²yinit)
    σ²x = isnothing(σ²xinit) ? ones(p) : copy(σ²xinit)
    κv = isnothing(κvinit) ? ones(p) : copy(κvinit)
    δw = isnothing(δwinit) ? zeros(qw) : copy(δwinit)
    δv = isnothing(δvinit) ? [zeros(qv[j]) for j in 1:p] : copy(δvinit)

    # pre-computation
    BwtBw = [Bw[k]' * Bw[k] for k in 1:K]
    BvxtBvx = [Bvx[j]' * Bvx[j] for j in 1:p]
    BvytBvy = [Bvy[j,k]' * Bvy[j,k] for j in 1:p, k in 1:K]

    # using design matrix approach
    βf = [b; β]
    A = [myA(k, n[k], K, id) for k in 1:K]
    Vf = [[A[k] stack(Bvy[:, k] .* δv)] for k in 1:K]
    resid = [z[k] - Vf[k] * βf - Bw[k] * δw for k in 1:K]

    # truncated limits
    zlower = [ifelse.(y[k] .== 0, nothing, 0) for k in 1:K]
    zupper = [ifelse.(y[k] .== 1, nothing, 0) for k in 1:K]

    # number of iterations
    niter = mcmcniter(nsamples, burnin, thin)
    nsamplesz = mcmcnsamples(niter, burnin, thinz)

    # initialize samples
    δw_samples = zeros(qw, nsamples + 1)
    δv_samples = [zeros(qv[i], nsamples + 1) for i in 1:p]
    z_samples = [zeros(n[k], nsamplesz + 1) for k in 1:K]
    βf_samples = zeros(pf, nsamples + 1)
    α_samples = zeros(p, nsamples + 1)
    σ²y_samples = zeros(K, nsamples + 1)
    σ²x_samples = zeros(p, nsamples + 1)
    κw_samples = zeros(1, nsamples + 1)
    κv_samples = zeros(p, nsamples + 1)

    δw_samples[:, 1] = δw
    [δv_samples[j][:, 1] = δv[j] for j in 1:p]
    [z_samples[k][:, 1] = z[k] for k in 1:K]
    βf_samples[:, 1] = βf
    σ²y_samples[:, 1] = σ²y
    σ²x_samples[:, 1] = σ²x
    κv_samples[:, 1] = κv

    # mcmc
    for i = 1:niter
        println("Iteration $i:")

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
            println("Iteration $i with saved sample $(saveid):")
            δw_samples[:, saveid] = δw
        end
        [resid[k] -= Bw[k] * δw for k in 1:K]

        # sample δv
        for j = 1:p
            Qv = sum(β[j]^2 ./ σ²y .* BvytBvy[j, :]) + BvxtBvx[j] / σ²x[j] + κv[j]*Pv[j]
            CQv = cholesky(Qv)
            [resid[k] += β[j] * Bvy[j, k] * δv[j] for k in 1:K]
            auxv = sum([σ²y[k]^(-1) * β[j] * Bvy[j,k]' * resid[k] for k in 1:K]) + σ²x[j]^(-1) * Bvx[j]' * (x[j])
            μv = CQv.UP \ (CQv.PtL \ auxv)
            δv[j] = μv + CQv.UP \ randn(qv[j])
            if saveiter
                δv_samples[j][:, saveid] = δv[j]
            end
            [resid[k] -= β[j] * Bvy[j, k] * δv[j] for k in 1:K]
        end

        # sample β
        Vf = [[A[k] stack(Bvy[:, k] .* δv)] for k in 1:K]
        Qᵦ = sum([σ²y[k]^(-1) * Vf[k]' * Vf[k] for k in 1:K]) + Σᵦ
        CQᵦ = cholesky(Qᵦ)
        μᵦ = CQᵦ.U \ (CQᵦ.L \ sum([σ²y[k]^(-1) * Vf[k]' * (z[k] - Bw[k] * δw) for k in 1:K]))
        βf = μᵦ + CQᵦ.U \ randn(pf)
        β = βf[K+1:end]
        if saveiter
            βf_samples[:, saveid] = βf
        end
        resid = [z[k] - Vf[k] * βf - Bw[k] * δw for k in 1:K]

        # sample σ²y
        a_σ²y = sigma2_a_prior .+ n / 2
        b_σ²y = [sigma2_b_prior + sum(resid[k].^2) / 2 for k in 1:K]
        σ²y = rand.(InverseGamma.(a_σ²y, b_σ²y))
        if saveiter
            σ²y_samples[:, saveid] = σ²y
        end

        # sample σ²x
        a_σ²x = sigma2_a_prior .+ m / 2
        b_σ²x = [sigma2_b_prior + sum((x[j] - Bvx[j]*δv[j]).^2) / 2 for j in 1:p]
        σ²x = rand.(InverseGamma.(a_σ²x, b_σ²x))
        if saveiter
            σ²x_samples[:, saveid] = σ²x
        end

        # sample κv
        a_κ = kappa_a_prior .+ qv / 2
        b_κ = kappa_b_prior .+ (transpose.(δv) .* Pv .* δv) ./ 2
        κv = rand.(Gamma.(a_κ, 1 ./ b_κ))
        if saveiter
            κv_samples[:, saveid] = κv
        end
    end

    Dict("w" => δw_samples, "v" => δv_samples, "beta" => βf_samples,
        "sigma2_y" => σ²y_samples, "sigma2_x" => σ²x_samples, "kappav" => κv_samples, "z" => z_samples)
end

