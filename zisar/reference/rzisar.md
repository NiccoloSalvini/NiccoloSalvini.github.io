# Simulate data from the ZISAR model

Draws one data set from the ZISAR model, on a regular lattice with queen
contiguity (the Monte Carlo design of the ZISAR paper) or on a weights
matrix supplied by the user. Covariates are independent standard normal;
the first element of `beta` and of `gamma` is the intercept. With
`sigma0 > 0` the inflation predictor has a SAR field of its own (the
general model of
[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md));
the default `sigma0 = 0` gives the restricted model.

## Usage

``` r
rzisar(
  n_side,
  beta,
  rho,
  gamma,
  kappa = 2,
  sigma = 0.5,
  rho0 = 0,
  sigma0 = 0,
  W = NULL
)
```

## Arguments

- n_side:

  Side of the square lattice, so \\n\\ = `n_side^2` units. Ignored when
  `W` is supplied.

- beta:

  Count-process coefficients, intercept first.

- rho:

  Spatial autoregressive parameter, in \\(-1, 1)\\.

- gamma:

  Zero-inflation coefficients, intercept first.

- kappa:

  Negative binomial shape.

- sigma:

  Standard deviation of the random effects.

- rho0, sigma0:

  Spatial autoregressive parameter and standard deviation of the
  inflation field (general model). With `sigma0 = 0` no field is drawn
  and the random-number stream is that of the restricted model.

- W:

  Optional row-standardised weights matrix (or `spdep` `listw`) to use
  instead of the lattice.

## Value

A list with

- data:

  Data frame with the counts `y`, the count covariates `x1, x2, ...` and
  the inflation covariates `z1, z2, ...`, ready for
  `zisar(y ~ x1 + ... | z1 + ..., data = data, W = W)`.

- W:

  The weights matrix (sparse).

- coords:

  Lattice coordinates (`row`, `col`), or `NULL` with a user-supplied
  `W`.

- truth:

  The parameters used.

- eta, g, mu, pi, s:

  Latent field, inflation field (zero in the restricted model), negative
  binomial means, inflation probabilities and structural-zero
  indicators.

## Examples

``` r
set.seed(42)
sim <- rzisar(15, beta = c(1, 0.5, -0.3), rho = 0.5, gamma = c(-1, 0.7))
head(sim$data)
#>   y      x1      x2      z1
#> 1 6  1.3710 -0.1740 -1.3038
#> 2 0 -0.5647  0.5157 -0.2509
#> 3 0  0.3631 -0.2344  0.1710
#> 4 0  0.6329 -0.6585 -0.4035
#> 5 7  0.4043  1.2502  0.1047
#> 6 4 -0.1061 -0.2718 -0.3189
mean(sim$data$y == 0)
#> [1] 0.4622
```
