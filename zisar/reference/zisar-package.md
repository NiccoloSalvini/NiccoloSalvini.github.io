# zisar: Zero-Inflated Spatial Autoregressive Models for Count Data

Fits the zero-inflated spatial autoregressive (ZISAR) model for spatial
count data with excess zeros. A simultaneous autoregressive (SAR)
process on the latent log-mean of a zero-inflated negative binomial
model carries spatial spillovers, and unit-level Gaussian random effects
absorb residual heterogeneity. An optional second SAR field in the
zero-inflation predictor captures spatially clustered structural zeros.
Parameters are estimated by maximising the Laplace approximation to the
marginal likelihood with 'TMB', whose automatic differentiation and
sparse matrix algebra make tens of thousands of spatial units feasible.
Methods are provided for direct, indirect and total marginal effects,
counterfactual prediction, residuals and simulation from the model.

## Details

The main function is
[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md),
which fits the model from a two-part formula
`y ~ count terms | zero-inflation terms` and a spatial weights matrix.
[`zisar_fit()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
is the matrix interface underneath it. Fitted models have the usual
methods
([summary()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[coef()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[vcov()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[logLik()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[predict()](https://niccolosalvini.github.io/zisar/reference/predict.zisar.md),
[residuals()](https://niccolosalvini.github.io/zisar/reference/residuals.zisar.md));
[`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md)
decomposes covariate effects into direct, indirect and total, and
[`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
simulates data from the model.

## References

Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
spatial autoregressive model for count data with random effects. Working
paper, Università Cattolica del Sacro Cuore.

## See also

Useful links:

- <https://github.com/NiccoloSalvini/zisar>

- Report bugs at <https://github.com/NiccoloSalvini/zisar/issues>

## Author

**Maintainer**: Niccolò Salvini <niccolo@binoocle.net>

Authors:

- Niccolò Salvini <niccolo@binoocle.net>

- Giuseppe Arbia

- Vincenzo Nardelli
