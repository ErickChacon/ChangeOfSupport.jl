# ChangeOfSupport.jl
"""
    _pseudointegral(x::Vector, b::RegularBsplines)

Return a dense matrix of the integral of the basis functions of `b` up to `x`.

It does not compute the integral of all the basis functions of `b`. It only computes those
the basis function numbers s-k+1:s where s is the last non-zero function of x.
"""

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

### Todo
