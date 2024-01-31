
# mcmcniter(nsamples::Int, burnin::Int, thin::Int) = burnin + 1 + thin * (nsamples - 1) + 1
# mcmcsave(i::Int, burnin::Int, thin::Int) = ((i - 1) > burnin) && (((i - 1) - burnin - 1) % thin == 0)
# mcmcid(i::Int, burnin::Int, thin::Int) = div((i - 1) - burnin - 1, thin) + 1 + 1

mcmcniter(nsamples::Int, burnin::Int, thin::Int) = burnin + 1 + thin * (nsamples - 1)
mcmcnsamples(niter::Int, burnin::Int, thin::Int) = length(burnin+1:thin:niter)
mcmcsave(i::Int, burnin::Int, thin::Int) = (i > burnin) && ((i - burnin - 1) % thin == 0)
mcmcid(i::Int, burnin::Int, thin::Int) = div(i - burnin - 1, thin) + 1

