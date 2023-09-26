# ChangeOfSupport.jl

Spatial modelling for data with different supports.

## Structure

### Composite Types

- [X] RegularKnots
- [X] RectilinearGrid
- [X] GMRF
    - [X] RGMRF
    - [X] CGMRF
- [X] RegularBsplines
- [X] NRegularBsplines

### Bayesian Inference

- [ ] Regular grid
- [ ] Rectilinear grid

### Todo

- [ ] Gibbs sampling
    - [ ] Take advantage of sparsity of the posterior precision matrix
- [ ] Integration with Turing.jl
- [ ] Compare the lines of cholesky when dense and when sparse:
```julia
# Qᵦ = Bw'Bw / 0.1^2 + P
# c1 = cholesky(Qᵦ)
# c2 = cholesky(sparse(Qᵦ))
```


### Todo compatibility with other packages

- [ ] Replace `RectilinearGrid` with the one in `Meshes.jl`.
- [ ] Need to move simulations to paper repository.


### PR

- [X] rangebars in `Makie.jl` should allow cycle colors.
- [ ] series in `Makie.jl` should not have labels.


