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

### Todo compatibility with other packages

- [ ] Replace `RectilinearGrid` with the one in `Meshes.jl`.
