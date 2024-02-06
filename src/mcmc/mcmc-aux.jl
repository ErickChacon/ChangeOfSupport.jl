# Function to create design matrix of intercept and biases for response source i with n
# observations and K response sources. We assume that the source id does not contain
# biases to have an identifiable model.
function myA(i, n, K, id)
    A = [ones(n) zeros(n, K-1)]
    if i != id
        A[:, 1+i-(i>id)] .= 1
    end
    A
end

# Functions to save samples depending of number of iterations, burninand thinning.
mcmcniter(nsamples::Int, burnin::Int, thin::Int) = burnin + 1 + thin * (nsamples - 1)
mcmcnsamples(niter::Int, burnin::Int, thin::Int) = length(burnin+1:thin:niter)
mcmcsave(i::Int, burnin::Int, thin::Int) = (i > burnin) && ((i - burnin - 1) % thin == 0)
mcmcid(i::Int, burnin::Int, thin::Int) = div(i - burnin - 1, thin) + 1

# Function to create initial values from MCMC chains
function getstate(mcmc)
    # dimensions
    K = length(mcmc["z"])

    # get current state
    init_b = mcmc["beta"][:, end][1:K]
    init_β = mcmc["beta"][:, end][K+1:end]
    init_σ²y = mcmc["sigma2_y"][:, end]
    init_σ²x = mcmc["sigma2_x"][:, end]
    init_κv = mcmc["kappav"][:, end]
    init_δw = mcmc["w"][:, end]
    init_δv = [x[:, end] for x in mcmc["v"]]

    # return tuple
    (binit = init_b, βinit = init_β, σ²yinit = init_σ²y, σ²xinit = init_σ²x,
        κvinit = init_κv, δwinit = init_δw, δvinit = init_δv)
end
