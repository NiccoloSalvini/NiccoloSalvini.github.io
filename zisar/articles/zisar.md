# Introduction to zisar

Spatial count data often have more zeros than a count distribution
allows, more variance than the mean, and dependence between neighbouring
units: the number of solar panels installed in a neighbourhood, of
emergency visits from a municipality, of firms of a sector in a
district. The ZISAR model handles the three at once, and **zisar**
estimates it by maximum likelihood.

``` r

library(zisar)
```

## The model

For $`n`$ spatial units linked by a weights matrix $`W`$
(row-standardised, with zero diagonal), the log-mean of the counts
follows a simultaneous autoregressive (SAR) process with unit-level
random effects,

``` math
\eta = \rho W \eta + X\beta + u, \qquad u \sim N(0, \sigma^2 I),
\qquad\text{so}\qquad \eta = (I - \rho W)^{-1}(X\beta + u).
```

Each count is a *structural zero* with probability $`\pi_i`$, and
otherwise negative binomial with mean $`\mu_i = \exp(\eta_i)`$ and shape
$`\kappa`$:

``` math
y_i \sim \pi_i \, 1\{y_i = 0\} + (1 - \pi_i)\,\mathrm{NB}(\mu_i, \kappa),
\qquad \mathrm{logit}(\pi_i) = z_i'\gamma.
```

The spatial multiplier $`(I - \rho W)^{-1}`$ acts on the *latent*
predictor, not on the observed counts, so the model has no endogenous
regressor and is estimated by maximum likelihood. Overdispersion has
three sources: the negative binomial shape $`\kappa`$, the structural
zeros $`\pi_i`$ and the random effects $`\sigma`$. The latent field is
integrated out with the Laplace approximation in
[TMB](https://cran.r-project.org/package=TMB).

This is the *restricted* model: the structural zeros depend on
covariates only. The *general* model adds a spatial field to the
inflation predictor, for zeros that cluster in space; see
[below](#clustered-structural-zeros).

## Fitting a model

[`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
simulates data on a regular lattice with queen contiguity. Here 400
units, a spatial parameter of 0.5, two covariates in the count process
and one in the zero-inflation process:

``` r

set.seed(2026)
sim <- rzisar(20, beta = c(0.5, 0.6, -0.4), rho = 0.5, gamma = c(-1, 0.9),
              kappa = 2, sigma = 0.4)
head(sim$data)
#>   y       x1      x2      z1
#> 1 1  0.52059  1.3369 -1.6053
#> 2 2 -1.07969  0.5527 -0.5025
#> 3 2  0.13924 -1.0565 -1.1959
#> 4 0 -0.08475 -0.1902 -1.0711
#> 5 6 -0.66664 -0.9092  0.7462
#> 6 0 -2.51609 -2.0009 -0.4369
mean(sim$data$y == 0)
#> [1] 0.465
```

[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
takes a two-part formula: the terms before `|` enter the count process,
those after it the zero-inflation process. `W` can be a base matrix, a
sparse `Matrix`, or an `spdep` `listw` object, with rows in the order of
the data.

``` r

fit <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W)
summary(fit)
#> 
#> Call:
#> zisar(formula = y ~ x1 + x2 | z1, data = sim$data, W = sim$W)
#> 
#> 400 spatial units, 46.5% zeros
#> 
#> Count model coefficients (SAR negative binomial, log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.3558     0.0905    3.93  8.5e-05 ***
#> x1            0.5918     0.0736    8.04  9.2e-16 ***
#> x2           -0.5495     0.0746   -7.37  1.7e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial, logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   -1.179      0.282   -4.18  2.9e-05 ***
#> z1             0.995      0.229    4.35  1.4e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Spatial and dispersion parameters:
#>       Estimate Std. Error z value Pr(>|z|)    
#> rho     0.6060     0.0972    6.24  4.5e-10 ***
#> kappa   2.0203     0.7320                     
#> sigma   0.3486     0.1223                     
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood (Laplace): -736.74 on 8 Df
#> AIC: 1489.5   BIC: 1521.4
#> Log-determinant: eigen.  Optimiser: converged after 30 iterations.
```

The estimates are close to the values used in the simulation
($`\beta = (0.5, 0.6, -0.4)`$, $`\gamma = (-1, 0.9)`$, $`\rho = 0.5`$,
$`\kappa = 2`$, $`\sigma = 0.4`$). The usual accessors work:

``` r

coef(fit)
#> count_(Intercept)          count_x1          count_x2  zero_(Intercept) 
#>            0.3558            0.5918           -0.5495           -1.1794 
#>           zero_z1               rho 
#>            0.9947            0.6060
sqrt(diag(vcov(fit)))
#> count_(Intercept)          count_x1          count_x2  zero_(Intercept) 
#>           0.09055           0.07363           0.07458           0.28232 
#>           zero_z1               rho 
#>           0.22871           0.09717
logLik(fit)
#> 'log Lik.' -736.7 (df=8)
AIC(fit)
#> [1] 1489
```

A formula without `|` uses the same terms in both processes, `y ~ x | 1`
gives a constant inflation probability, and an
[`offset()`](https://rdrr.io/r/stats/offset.html) term in the count part
enters the log-mean with coefficient one, for exposures:
`y ~ x1 + offset(log(households)) | z1`.

## Spatial weights from real geography

With polygons, `spdep` builds the weights; pass the `listw` object (or
the matrix) to
[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md).
Every unit must have a row in the data: units cannot be dropped, because
they all enter the spatial multiplier.

``` r

library(sf)
library(spdep)
areas <- st_read("areas.gpkg")
nb <- poly2nb(areas, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
fit <- zisar(y ~ x1 + x2 | z1, data = areas, W = lw)
```

## Direct, indirect and total effects

In a spatial model a coefficient is not the effect of its covariate: a
change in one unit moves the expected count there (direct effect) and,
through the spatial multiplier, in the other units (indirect effect, or
spillover).
[`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md)
averages them over units, with simulation standard errors:

``` r

zisar_effects(fit, B = 100)
#> ZISAR average effects on E[y] (simulation SEs, B = 100)
#> 
#>              Direct         Indirect            Total
#> x1  1.6312 (0.2927)  2.2126 (0.9602)  3.8438 (1.1102)
#> x2 -1.5147 (0.2520) -2.0545 (0.7905) -3.5692 (0.8750)
#> z1 -0.5911 (0.1170)  0.0000 (0.0000) -0.5911 (0.1170)
```

`z1` acts only through the probability of a structural zero, which
concerns the unit itself: its indirect effect is zero. A covariate
entering both processes (the same name in both parts of the formula)
combines the two channels.

## Counterfactual scenarios

[`predict()`](https://rdrr.io/r/stats/predict.html) gives the fitted
means (`type = "response"`), the means of the count component
(`"count"`), the inflation probabilities (`"zero"`) or the latent field
(`"eta"`). With `newdata` for the same units it computes a scenario:
here `x1` rises by one in a single unit, and the expected counts of its
neighbours rise as well.

``` r

nd <- sim$data
nd$x1[190] <- nd$x1[190] + 1
change <- predict(fit, newdata = nd) - predict(fit)
round(change[c(190, which(sim$W[190, ] > 0))], 3)
#> [1] 1.615 0.064 0.083 0.016 0.193 0.211 0.256 0.264 0.054
sum(change)
#> [1] 3.594
```

## Comparing models

`fixed` holds parameters at given values. With $`\rho = 0`$ and
$`\sigma`$ close to zero the model reduces to the aspatial zero-inflated
negative binomial, the same model as `pscl::zeroinfl(dist = "negbin")`:

``` r

fit_zinb <- zisar(y ~ x1 + x2 | z1, data = sim$data, W = sim$W,
                  fixed = list(rho = 0, sigma = 0.001))
AIC(fit_zinb, fit)
#>          df  AIC
#> fit_zinb  6 1507
#> fit       8 1489
```

The test of $`\rho = 0`$ is in the summary (Wald). The aspatial model
with random effects (`fixed = list(rho = 0)`) is nested in ZISAR, but on
spatially dependent data it is often not estimable: the random effects
have to absorb the spatial structure, $`\sigma`$ grows, and the Laplace
approximation breaks down (see below).

## Clustered structural zeros

In the restricted model, structural zeros depend on $`z_i`$ only. When
whole areas are out of the market (a planning restriction, no supplier
nearby), the zeros cluster beyond what the covariates explain. The
general model gives the inflation predictor a SAR field of its own,

``` math
\mathrm{logit}(\pi_i) = z_i'\gamma + g_i, \qquad
g = (I - \rho_0 W)^{-1} e_0, \qquad e_0 \sim N(0, \sigma_0^2 I),
```

independent of $`u`$.
[`rzisar()`](https://niccolosalvini.github.io/zisar/reference/rzisar.md)
draws it when `sigma0 > 0`, and `zero_spatial = TRUE` fits it:

``` r

set.seed(2)
sim0 <- rzisar(15, beta = c(1, 0.5), rho = 0.4, gamma = c(-1, 0.8),
               rho0 = 0.6, sigma0 = 1)
fit_r <- zisar(y ~ x1 | z1, data = sim0$data, W = sim0$W)
fit_g <- zisar(y ~ x1 | z1, data = sim0$data, W = sim0$W,
               zero_spatial = TRUE)
summary(fit_g)$spatial
#>        Estimate Std. Error z value  Pr(>|z|)
#> rho      0.4349     0.2585   1.683 9.246e-02
#> kappa    1.7261     0.3713      NA        NA
#> sigma    0.0010     0.3483      NA        NA
#> rho0     0.8538     0.1289   6.624 3.486e-11
#> sigma0   0.4795     0.2994      NA        NA
```

The restricted model is the general one with $`\sigma_0 = 0`$, on the
boundary of the parameter space, so the likelihood-ratio statistic has a
mixture distribution and the $`\chi^2_2`$ p-value is conservative:

``` r

lr <- 2 * (logLik(fit_g) - logLik(fit_r))
pchisq(as.numeric(lr), df = 2, lower.tail = FALSE)
#> [1] 0.1966
```

The two latent fields compete for the same zeros (a zero is either a
structural zero, with a large $`g_i`$, or a negative binomial zero, with
a small $`\eta_i`$), so the general model is harder to estimate and more
exposed to the breakdown described next. Fit the restricted model first,
and read $`\rho_0`$ and $`\sigma_0`$ only from a general fit with
`laplace_ok = TRUE`.

## When the Laplace approximation breaks down

For a zero count, the conditional distribution of $`\eta_i`$ can be
bimodal: the zero is either structural (and $`\mu_i`$ can be large) or
negative binomial (and $`\mu_i`$ is small). When $`\sigma`$ is large
this produces spurious optima, with a huge gradient and sometimes a
positive log-likelihood. It happens mostly in small samples and with
misspecified models, for example with an omitted covariate.
[`zisar()`](https://niccolosalvini.github.io/zisar/reference/zisar.md)
checks every fit, restarts from other values when a check fails, and
warns if no usable optimum is found; `fit$laplace_ok` records the
outcome and [`summary()`](https://rdrr.io/r/base/summary.html) reports
it.

``` r

fit$laplace_ok
#> [1] TRUE
fit$max_grad
#> [1] 0.000214
```

## Large data sets

The log-determinant $`\log|I - \rho W|`$ is computed from the
eigenvalues of $`W`$ for up to 5000 units, and from a sparse Cholesky
factorisation beyond (`logdet = "sparse"`), which keeps memory and time
manageable for tens of thousands of units.
[`zisar_effects()`](https://niccolosalvini.github.io/zisar/reference/zisar_effects.md)
never forms the dense $`n \times n`$ multiplier either.

``` r

fit <- zisar(y ~ x1 + x2 | z1, data = big, W = W_big, logdet = "sparse")
```

## Reference

Salvini, N., Arbia, G. and Nardelli, V. (2026). ZISAR: A zero-inflated
spatial autoregressive model for count data with random effects. Working
paper, Università Cattolica del Sacro Cuore.
