# Fit a zero-inflated spatial autoregressive (ZISAR) model

Fits the ZISAR model of Salvini, Arbia and Nardelli (2026): a
zero-inflated negative binomial model whose log-mean follows a
simultaneous autoregressive (SAR) process with unit-level random
effects. By default the inflation probability depends on covariates only
(the restricted model); with `zero_spatial = TRUE` it has a spatial
random field of its own (the general model). `zisar()` takes a two-part
formula and a data frame; `zisar_fit()` is the matrix interface
underneath it.

## Usage

``` r
zisar(
  formula,
  data,
  W,
  zero_spatial = FALSE,
  fixed = NULL,
  start = NULL,
  logdet = c("auto", "eigen", "sparse"),
  control = list(),
  silent = TRUE
)

zisar_fit(
  y,
  X,
  Z = NULL,
  W,
  offset = NULL,
  zero_spatial = FALSE,
  fixed = NULL,
  start = NULL,
  logdet = c("auto", "eigen", "sparse"),
  control = list(),
  silent = TRUE
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: the terms before `|` enter
  the count process, those after it the zero-inflation process. Without
  `|`, both processes use the same terms. Use `| 1` for a constant
  inflation probability. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the count
  part enters the log-mean with coefficient one.

- data:

  A data frame with one row per spatial unit, in the order of the rows
  of `W`. Missing values are not allowed: every unit enters the spatial
  multiplier, so none can be dropped.

- W:

  Spatial weights: an \\n \times n\\ matrix (base or sparse 'Matrix') or
  an `spdep` `listw` object, with zero diagonal and preferably
  row-standardised.

- zero_spatial:

  If `TRUE`, fit the general model, with a SAR random field in the
  inflation predictor (see Model).

- fixed:

  Optional named list of parameters held at fixed values: `rho`,
  `kappa`, `sigma`, and in the general model `rho0` and `sigma0`
  (scalars, natural scale), `beta` (length `ncol(X)`) or `gamma` (length
  `ncol(Z)`). In `beta` and `gamma`, `NA` entries are estimated.

- start:

  Optional named list of starting values, with the same names and scales
  as `fixed`. Defaults: `beta` from a Poisson regression, inflation
  intercept from the share of zeros, `rho = 0`, `kappa = 1`,
  `sigma = 0.5`, and `rho0 = 0`, `sigma0 = 0.5`.

- logdet:

  Method for the log-determinant: `"auto"`, `"eigen"` or `"sparse"`. See
  Estimation.

- control:

  List of control parameters for
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), merged with
  the defaults `eval.max = 1000`, `iter.max = 500`, `rel.tol = 1e-10`.

- silent:

  If `FALSE`, print the 'TMB' tracing output.

- y:

  Vector of non-negative integer counts.

- X:

  Design matrix of the count process, \\n \times p\\; include a column
  of ones for an intercept.

- Z:

  Design matrix of the zero-inflation process, \\n \times q\\. `NULL`
  gives an intercept only.

- offset:

  Optional vector of length \\n\\ added to the count log-mean (for
  example the log of an exposure).

## Value

An object of class `"zisar"`, a list with, among others:

- coefficients:

  List with the `count` (\\\beta\\) and `zero` (\\\gamma\\)
  coefficients.

- rho, kappa, sigma:

  Spatial autoregressive parameter, negative binomial shape and
  random-effect standard deviation.

- rho0, sigma0:

  In the general model, the spatial autoregressive parameter and
  standard deviation of the inflation field; `NULL` otherwise.

- eta:

  Mode of the latent field \\\eta\\ given the data.

- g:

  Mode of the inflation field \\g\\ (zero in the restricted model).

- fitted.values:

  Fitted means \\E\[y_i\] = (1 - \pi_i)\mu_i\\.

- loglik:

  Laplace-approximated marginal log-likelihood.

- convergence, message:

  Convergence code and message of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- max_grad, laplace_ok:

  Largest gradient component at the optimum, and whether the optimum
  passed the checks described in Estimation.

- tmb:

  The 'TMB' object, for advanced use.

Use the methods
[summary()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[coef()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[vcov()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[logLik()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md),
[predict()](https://niccolosalvini.github.io/zisar/reference/predict.zisar.md),
[residuals()](https://niccolosalvini.github.io/zisar/reference/residuals.zisar.md)
and
[`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md)
rather than the components directly.

## Model

For units \\i = 1, \dots, n\\ linked by a spatial weights matrix \\W\\,
\$\$\eta = \rho W \eta + X\beta + u, \qquad u \sim N(0, \sigma^2 I),\$\$
so that \\\eta = (I - \rho W)^{-1}(X\beta + u)\\. Each count is zero
with probability \\\pi_i\\ (a structural zero) and otherwise negative
binomial, \$\$y_i \sim \pi_i \\ 1\\y_i = 0\\ + (1 - \pi_i)\\
\mathrm{NB}(\mu_i, \kappa), \qquad \mu_i = \exp(\eta_i + o_i), \qquad
\mathrm{logit}(\pi_i) = z_i'\gamma,\$\$ with \\\mathrm{Var}(y_i \mid
\text{not a structural zero}) = \mu_i + \mu_i^2/\kappa\\ and \\o_i\\ an
optional offset (log exposure). The offset scales unit \\i\\ only: it
does not pass through the spatial multiplier.

In the general model (`zero_spatial = TRUE`) the inflation predictor has
a SAR random field of its own, \$\$\mathrm{logit}(\pi_i) = z_i'\gamma +
g_i, \qquad g = \rho_0 W g + e_0, \qquad e_0 \sim N(0, \sigma_0^2
I),\$\$ independent of \\u\\, so that structural zeros can cluster in
space beyond what the covariates explain. The restricted model is the
limit \\\sigma_0 \to 0\\.

## Estimation

The latent field \\\eta\\ is integrated out with the Laplace
approximation in 'TMB' (together with \\g\\ in the general model), and
the marginal likelihood is maximised with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) over \\(\beta,
\gamma, \rho, \kappa, \sigma)\\, plus \\(\rho_0, \sigma_0)\\ in the
general model, with \\\rho = \tanh(\rho\_{raw})\\ restricted to \\(-1,
1)\\. Standard errors come from the inverse Hessian of the Laplace
objective; those of \\\rho\\, \\\kappa\\, \\\sigma\\, \\\rho_0\\ and
\\\sigma_0\\ from the delta method.

The Laplace approximation can break down. For a zero count, the
conditional distribution of \\\eta_i\\ can be bimodal (a structural zero
with a large mean, or a negative binomial zero with a small one), and
when \\\sigma\\ is large the approximation produces spurious optima with
a positive log-likelihood and a huge gradient. This happens mostly with
small samples and misspecified models (an omitted covariate inflates
\\\sigma\\). Each fit is therefore checked: when the optimiser does not
converge, or the optimum has a gradient above 0.01 or a positive
log-likelihood, the fit is restarted from three other values of
\\\kappa\\ and \\\sigma\\ and from up to three random perturbations of
the best point (this uses the random number generator). If no usable
optimum is found, `zisar()` warns and `laplace_ok` is `FALSE`.

The general model is more exposed to this breakdown: a zero count can
come from a large \\g_i\\ (a structural zero) or from a small \\\eta_i\\
(a negative binomial zero), and the two latent fields compete for the
same zeros. Fit the restricted model first, and treat a general fit with
`laplace_ok = FALSE` as uninformative about \\(\rho_0, \sigma_0)\\.

The log-determinant \\\log\|I - \rho W\|\\ is computed either from the
eigenvalues of \\W\\, found once in \\O(n^3)\\ time and \\O(n^2)\\
memory (`logdet = "eigen"`), or from a sparse Cholesky factorisation of
\\(I - \rho W)'(I - \rho W)\\ at each evaluation (`logdet = "sparse"`),
which scales to tens of thousands of units. `"auto"` uses the sparse
method when \\n \> 5000\\.

## Nested models

`fixed` holds parameters at given values, which gives the special cases
of the model: `fixed = list(rho = 0)` is the aspatial zero-inflated
negative binomial with random effects (ZINB-RE), and with an
intercept-only inflation part (`y ~ x | 1`), `fixed = list(gamma = -20)`
sets \\\pi_i \approx 0\\, the SAR negative binomial with random effects
(SAR-NB-RE). Likelihood-ratio tests between nested fits use
[logLik()](https://niccolosalvini.github.io/zisar/reference/zisar-methods.md).
The restricted model sits on the boundary of the general one (\\\sigma_0
= 0\\), so the likelihood-ratio statistic between them does not have the
usual chi-squared distribution: its p-value from \\\chi^2_2\\ is
conservative.

## References

Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
spatial autoregressive model for count data with random effects. Working
paper, Università Cattolica del Sacro Cuore.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1–21.
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05)

## See also

[`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md)
for direct, indirect and total effects,
[`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
to simulate from the model.

## Examples

``` r
set.seed(1)
sim <- rzisar(12, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8))
fit <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W)
summary(fit)
#> 
#> Call:
#> zisar(formula = y ~ x1 | z1, data = sim$data, W = sim$W)
#> 
#> 144 spatial units, 38.9% zeros
#> 
#> Count model coefficients (SAR negative binomial, log link):
#>             Estimate Std. Error z value     Pr(>|z|)    
#> (Intercept)    1.622      0.809    2.01        0.045 *  
#> x1             0.669      0.114    5.87 0.0000000045 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial, logit link):
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)   -0.805      0.264   -3.05   0.0023 **
#> z1             0.652      0.224    2.91   0.0037 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Spatial and dispersion parameters:
#>       Estimate Std. Error z value Pr(>|z|)
#> rho      0.106      0.412    0.26      0.8
#> kappa    1.461      0.803                 
#> sigma    0.203      0.628                 
#> 
#> Log-likelihood (Laplace): -344.3 on 7 Df
#> AIC: 702.61   BIC: 723.39
#> Log-determinant: eigen.  Optimiser: converged after 24 iterations.
#> 

# Aspatial ZINB with random effects, and a likelihood-ratio test of rho = 0
fit0 <- zisar(y ~ x1 | z1, data = sim$data, W = sim$W,
              fixed = list(rho = 0))
lr <- 2 * (logLik(fit) - logLik(fit0))
pchisq(as.numeric(lr), df = 1, lower.tail = FALSE)
#> [1] 0.8091

# \donttest{
# General model: zeros clustered by a spatial field of their own
set.seed(2)
sim2 <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
               rho0 = 0.6, sigma0 = 1)
fit2 <- zisar(y ~ x1 | z1, data = sim2$data, W = sim2$W,
              zero_spatial = TRUE)
summary(fit2)
#> 
#> Call:
#> zisar(formula = y ~ x1 | z1, data = sim2$data, W = sim2$W, zero_spatial = TRUE)
#> 
#> 225 spatial units, 43.1% zeros
#> 
#> Count model coefficients (SAR negative binomial, log link):
#>             Estimate Std. Error z value   Pr(>|z|)    
#> (Intercept)   0.9574     0.4441    2.16      0.031 *  
#> x1            0.4343     0.0836    5.19 0.00000021 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial, logit link):
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)   -0.765      0.317   -2.41    0.016 * 
#> z1             0.717      0.232    3.09    0.002 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Spatial and dispersion parameters (rho0, sigma0: inflation field):
#>        Estimate Std. Error z value       Pr(>|z|)    
#> rho       0.435      0.258    1.68          0.092 .  
#> kappa     1.726      0.371                           
#> sigma     0.001      0.348                           
#> rho0      0.854      0.129    6.62 0.000000000035 ***
#> sigma0    0.480      0.299                           
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood (Laplace): -499.55 on 9 Df
#> AIC: 1017.1   BIC: 1047.8
#> Log-determinant: eigen.  Optimiser: converged after 41 iterations.
#> 
# }
```
