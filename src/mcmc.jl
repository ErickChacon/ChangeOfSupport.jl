
# Mcmc sampling

for filename in ["mcmc-aux.jl", "mcmc.jl", "mcmc-model.jl", "mcmc-model-y.jl",
    "mcmc-model-standard.jl", "mcmc-model-standard-no-w.jl"]
    include(joinpath("mcmc", filename))
end

