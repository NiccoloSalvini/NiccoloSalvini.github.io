# Changelog

## zisar 0.1.0

First release.

- [`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
  fits the ZISAR model from a two-part formula
  `y ~ count terms | zero-inflation terms` and a spatial weights matrix
  (base matrix, sparse `Matrix` or `spdep` `listw`);
  [`zisar_fit()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
  is the matrix interface.
- Estimation by Laplace approximation in ‘TMB’, with the log-determinant
  from the eigenvalues of `W` or, for large `n`, from a sparse Cholesky
  factor.
- `fixed` holds parameters at given values, for the nested ZINB-RE
  (`rho = 0`) and SAR-NB-RE models.
- `zero_spatial = TRUE` fits the general model, with a SAR random field
  in the zero-inflation predictor (parameters `rho0` and `sigma0`);
  [`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
  simulates from it with `sigma0 > 0`.
- [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the count
  part, for exposures.
- Methods: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (hence
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html)),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) (including
  counterfactual scenarios on the same units) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).
- [`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md):
  average direct, indirect and total effects with simulation standard
  errors, including the zero-inflation channel.
- [`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
  simulates from the model;
  [`lattice_weights()`](https://niccolosalvini.github.io/zisar/reference/lattice_weights.md)
  builds lattice contiguity weights.
