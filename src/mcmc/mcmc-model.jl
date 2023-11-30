function myA(i, n, K, id)
    A = [ones(n) zeros(n, K-1)]
    if i != id
        A[:, 1+i-(i>id)] .= 1
    end
    A
end

function sample_model(y, x, Bw, Bvx, Bvy, Pw, Pv, σ²y, σ²x, κw, κv, β₀, β, b, α, δw, δv, id; niter = 10)

    α = copy(α)

    # dimensions
    K = length(y)
    p = length(x)
    pf = K + p
    n = length.(y)
    m = map(x -> size(x, 1), x)
    qw = size(Bw[1], 2)
    qv = map(x -> size(x, 2), Bvx)

    # # hyperparameters
    # sigma2_a_prior = 0.001
    # sigma2_b_prior = 0.001
    # kappa_a_prior = 0.001
    # kappa_b_prior = 0.001
    Σᵦ = 3I
    σ²α = repeat([1], p)

    # initialize samples
    δw_samples = zeros(qw, niter)
    δv_samples = [zeros(qv[i], niter) for i in 1:p]
    βf_samples = zeros(pf, niter)
    α_samples = zeros(p, niter)
    σ²y_samples = zeros(K, niter)
    σ²x_samples = zeros(p, niter)
    κw_samples = zeros(1, niter)
    κv_samples = zeros(p, niter)

    # σ²_samples[1, 1] = σ²
    # κ_samples[1, 1] = κ

    # pre-computation
    BwtBw = [Bw[k]' * Bw[k] for k in 1:K]
    BvxtBvx = [Bvx[j]' * Bvx[j] for j in 1:p]
    BvytBvy = [Bvy[j,k]' * Bvy[j,k] for j in 1:p, k in 1:K]

    # using design matrix approach
    βf = [β₀; b; β]
    A = [myA(k, n[k], K, id) for k in 1:K]
    Vf = [[A[k] stack(Bvy[:, k] .* δv)] for k in 1:K]
    resid = [y[k] - Vf[k] * βf - Bw[k] * δw for k in 1:K]

    # # mcmc
    for i = 2:niter

        # sample δw
        Qw = sum(BwtBw ./ σ²y) + κw * Pw
        CQw = cholesky(Qw)
        [resid[k] += Bw[k] * δw for k in 1:K]
        aux = sum([σ²y[k]^(-1) * Bw[k]' * resid[k] for k in 1:K])
        μw =  CQw.UP \ (CQw.PtL \ aux)
        δw = μw + CQw.UP \ randn(qw)
        δw_samples[:, i] = δw
        [resid[k] -= Bw[k] * δw for k in 1:K]

        # sample δv
        for j = 1:p
            Qv = sum(β[j]^2 ./ σ²y .* BvytBvy[j, :]) + BvxtBvx[j] / σ²x[j] + κv[j]*Pv[j]
            CQv = cholesky(Qv)
            [resid[k] += β[j] * Bvy[j, k] * δv[j] for k in 1:K]
            auxv = sum([σ²y[k]^(-1) * β[j] * Bvy[j,k]' * resid[k] for k in 1:K]) + σ²x[j]^(-1) * Bvx[j]' * (x[j] .- α[j])
            μv = CQv.UP \ (CQv.PtL \ auxv)
            δv[j] = μv + CQv.UP \ randn(qv[j])
            δv_samples[j][:, i] = δv[j]
            [resid[k] -= β[j] * Bvy[j, k] * δv[j] for k in 1:K]
        end

        # sample β
        Vf = [[A[k] stack(Bvy[:, k] .* δv)] for k in 1:K]
        Qᵦ = sum([σ²y[k]^(-1) * Vf[k]' * Vf[k] for k in 1:K]) + Σᵦ
        CQᵦ = cholesky(Qᵦ)
        μᵦ = CQᵦ.U \ (CQᵦ.L \ sum([σ²y[k]^(-1) * Vf[k]' * (y[k] - Bw[k] * δw) for k in 1:K]))
        βf = μᵦ + CQᵦ.U \ randn(pf)
        β = βf[K+1:end]
        βf_samples[:, i] = βf
        resid = [y[k] - Vf[k] * βf - Bw[k] * δw for k in 1:K]

        # sample α
        for j = 1:p
            Qα = σ²x[j]^(-1) * m[j] + σ²α[j]
            μsα = (σ²x[j]^(-1) * sum(x[j] - Bvx[j]*δv[j]))
            α[j] = Qα^(-1) * (μsα + randn())
            α_samples[j, i] = α[j]
        end

        # # sample σ²
        # a_σ² = sigma2_a_prior + ny / 2
        # b_σ² = sigma2_b_prior + sum((y - Bw * β).^2) / 2
        # σ² = rand(InverseGamma(a_σ², b_σ²))
        # σ²_samples[1, i] = σ²
        #
        # if !fix_hyper
        #     # sample κ
        #     a_κ = kappa_a_prior + (p - 1) / 2
        #     b_κ = kappa_b_prior + transpose(β) * P * β / 2
        #     κ = rand(Gamma(a_κ, 1 / b_κ))
        # end
        # κ_samples[1, i] = κ
    end

    # δw_samples, σ²_samples, κ_samples
    # out
    δw_samples, δv_samples, βf_samples, α_samples
    # A, V, Vf
    # BvytBvy
end

