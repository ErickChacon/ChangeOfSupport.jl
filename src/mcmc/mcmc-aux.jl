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

