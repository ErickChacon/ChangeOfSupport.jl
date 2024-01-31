
# Mcmc sampling

for filename in ["mcmc-aux.jl", "mcmc.jl", "mcmc-model.jl", "mcmc-model-y.jl"]
    include(joinpath("mcmc", filename))
end

