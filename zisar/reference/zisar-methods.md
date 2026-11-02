# Methods for fitted ZISAR models

Standard accessors for objects returned by
[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
and
[`zisar_fit()`](https://niccolosalvini.github.io/zisar/reference/zisar.md).

## Usage

``` r
# S3 method for class 'zisar'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zisar'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zisar'
logLik(object, ...)

# S3 method for class 'zisar'
nobs(object, ...)

# S3 method for class 'zisar'
fitted(object, ...)

# S3 method for class 'zisar'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zisar'
summary(object, ...)

# S3 method for class 'summary.zisar'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object, x:

  A fitted `"zisar"` object.

- model:

  Which coefficients: `"full"` (count, zero-inflation, \\\rho\\ and, in
  the general model, \\\rho_0\\), `"count"` or `"zero"`.

- ...:

  Not used.

- digits:

  Number of significant digits to print.

- signif.stars:

  Logical: show significance stars? See
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix,
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"`, [`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
spatial units, [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
the fitted means \\E\[y_i\] = (1 - \pi_i)\mu_i\\,
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`"summary.zisar"`.

## Details

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) cover the regression
coefficients and \\\rho\\ (and \\\rho_0\\ in the general model); the
negative binomial shape \\\kappa\\ and the random-effect standard
deviations \\\sigma\\ and \\\sigma_0\\ are reported by
[`summary()`](https://rdrr.io/r/base/summary.html) and stored in
`object$kappa`, `object$sigma` and `object$sigma0`. The covariance of
\\\rho\\ and \\\rho_0\\ comes from the delta method; rows and columns of
parameters held fixed are `NA`.

[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the
Laplace-approximated marginal log-likelihood, with degrees of freedom
equal to the number of estimated parameters, so
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work on fitted
models.

## Examples

``` r
set.seed(1)
sim <- rzisar(10, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
coef(fit)
#> count_(Intercept)          count_x1  zero_(Intercept)           zero_z1 
#>            1.3107            0.3799           -1.0503            0.6185 
#>               rho 
#>            0.2311 
sqrt(diag(vcov(fit)))
#> count_(Intercept)          count_x1  zero_(Intercept)           zero_z1 
#>            0.8154            0.1575            0.4006            0.3068 
#>               rho 
#>            0.4686 
AIC(fit)
#> [1] 517.9
```
