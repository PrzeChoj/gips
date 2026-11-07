# The constructor of a `gips` class.

Create a `gips` object. This object will contain initial data and all
other information needed to find the most likely invariant permutation.
It will not perform optimization; one must call the
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
function to do it. See the examples below.

## Usage

``` r
gips(
  S,
  number_of_observations,
  delta = 3,
  D_matrix = NULL,
  was_mean_estimated = TRUE,
  perm = ""
)

new_gips(
  list_of_gips_perm,
  S,
  number_of_observations,
  delta,
  D_matrix,
  was_mean_estimated,
  optimization_info
)

validate_gips(g)
```

## Arguments

- S:

  A p by p empirical covariance matrix. This is the standard,
  single-sample use of `gips()`. For the optional multi-sample
  interface, `S` can also be a list of covariance matrices; see the
  **Multi-sample** section below. Row and column names are optional. If
  both are supplied, they should identify the same variables in the same
  order.

      When `Z` is the observed data:

  - if one does not know the theoretical mean and has to estimate it
    with the observed mean, use `S = cov(Z)`, and leave parameter
    `was_mean_estimated = TRUE` as default;

  - if one knows the theoretical mean is 0, use
    `S = crossprod(Z) / number_of_observations`, and set parameter
    `was_mean_estimated = FALSE`.

- number_of_observations:

  A number of data points that `S` is based on. For multi-sample input,
  use one value per matrix in `S`; see the **Multi-sample** section
  below.

- delta:

  A number, hyper-parameter of a Bayesian model. It has to be strictly
  bigger than 1. See the **Hyperparameters** section below. For
  multi-sample input, it can also be a numeric vector with one value per
  matrix in `S`; see the **Multi-sample** section below.

- D_matrix:

  Hyper-parameter of a Bayesian model. Either `NULL` (default), a
  symmetric, positive-definite matrix of the same size as `S` for
  single-sample input, or a list of such matrices (one per matrix in
  `S`) for multi-sample input. List elements cannot be `NULL`. When
  `NULL`, defaults are derived from the data separately for each sample.
  For more details, see the **Hyperparameters** and **Multi-sample**
  sections below.

- was_mean_estimated:

  A boolean.

  - Set `TRUE` (default) when your `S` parameter is a result of a
    [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function.

  - Set `FALSE` when your `S` parameter is a result of a
    `crossprod(Z) / number_of_observations` calculation.

- perm:

  An optional permutation to be the base for the `gips` object. It can
  be of a `gips_perm` or a `permutation` class, or anything the function
  [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html)
  can handle. It can also be of a `gips` class, but it will be
  interpreted as the underlying `gips_perm`.

- list_of_gips_perm:

  A list with a single element of a `gips_perm` class. The base object
  for the `gips` object.

- optimization_info:

  For internal use only. `NULL` or the list with information about the
  optimization process.

- g:

  Object to be checked whether it is a proper object of a `gips` class.

## Value

`gips()` returns an object of a `gips` class after the safety checks.

`new_gips()` returns an object of a `gips` class without safety checks.

`validate_gips()` returns its argument unchanged. If the argument is not
a proper element of a `gips` class, it produces an error.

## Functions

- `new_gips()`: Constructor. It is only intended for low-level use.

- `validate_gips()`: Validator. It is only intended for low-level use.

## Methods for a `gips` class

- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md)

- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md)

- [`print.gips()`](https://przechoj.github.io/gips/reference/print.gips.md)

- [`logLik.gips()`](https://przechoj.github.io/gips/reference/logLik.gips.md)

- [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)

- [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)

- [`as.character.gips()`](https://przechoj.github.io/gips/reference/as.character.gips.md)

## Hyperparameters

We encourage you to try `D_matrix = d * I`, where `I` is a `p`
\\\times\\ `p` identity matrix and `d > 0` for some different `d`.
Smaller `d` favor models with more equality constraints and fewer free
covariance parameters. This tendency depends on the data and `delta`.
For fixed dimension, sample size and `delta`, as `d` tends to infinity,
the posterior distribution over symmetry groups converges to a
distribution independent of `S`.

Using `S` with `D_matrix = d * I` gives the same results as using
`S / d` with `D_matrix = I`, with the other parameters unchanged.

The default for `D_matrix` is `D_matrix = d * I`, where
`d = mean(diag(S))`. For group selection, this is equivalent to using
`S / d` with `D_matrix = I`, so that the mean diagonal entry of `S / d`
is 1.

The Bayesian model uses a [Wishart-type
prior](https://en.wikipedia.org/wiki/Wishart_distribution) on the
precision matrix (the inverse covariance matrix).

The hyperparameters affect symmetry selection. For a fixed permutation,
[`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
does not use `delta` or `D_matrix`.

For a brief introduction, see the **Bayesian model selection** section
in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or in its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html).

For analysis of the Hyperparameters influence, see **Section 3.2.** of
"Learning Permutation Symmetry of a Gaussian Vector with gips in R" by
`gips` developers Adam Chojecki, Paweł Morgen, and Bartosz Kołodziejek,
Journal of Statistical Software;
[doi:10.18637/jss.v112.i07](https://doi.org/10.18637/jss.v112.i07) .

## Multi-sample

The usual and recommended starting point is the single-sample call
`gips(S, number_of_observations)`.

The optional multi-sample interface is for G independent groups that are
assumed to share the same permutation symmetry, while each group keeps
its own covariance matrix and sample size. Pass a list of G covariance
matrices as `S` and a numeric vector of length G as
`number_of_observations`:

    S_by_group <- list(control = S1, treatment = S2, follow_up = S3)
    g <- gips(S_by_group, c(n1, n2, n3))

All matrices in `S` must be square covariance matrices of the same size,
and their rows and columns must refer to the same variables in the same
order. `gips()` applies one shared permutation to matching matrix
indices; it does not match or reorder variables using matrix names.
Names of the list elements can identify groups, while matching row and
column names on the matrices can identify variables. These names are
used as labels in multi-sample heatmaps in
[`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md).
`D_matrix` should then be a list of G positive-definite matrices
(defaulting to `diag(mean(diag(S_g)), p)` for each group). `delta` can
be a scalar (broadcast to all groups) or a vector of length G (one value
per group); it defaults to `rep(3, G)`. `was_mean_estimated` is
intentionally one boolean value applied to all groups, not a per-group
setting.

## See also

- [`stats::cov()`](https://rdrr.io/r/stats/cor.html) - The `S`
  parameter, as an empirical covariance matrix, is most of the time a
  result of the [`cov()`](https://rdrr.io/r/stats/cor.html) function.
  For more information, see [Wikipedia - Estimation of covariance
  matrices](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices).

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The function that finds the Maximum A Posteriori (MAP) Estimator for a
  given `gips` object.

- [`gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md) -
  The constructor of a `gips_perm` class. The `gips_perm` object is used
  as the base object for the `gips` object. To be more precise, the base
  object for `gips` is a one-element list of a `gips_perm` object.

## Examples

``` r
require("MASS") # for mvrnorm()

perm_size <- 5
mu <- runif(5, -10, 10) # Assume we don't know the mean
sigma_matrix <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)

g_map <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
g_map
#> The permutation (1,2,3,4,5):
#>  - was found after 67 posterior calculations;
#>  - is 953.095 times more likely than the () permutation.

summary(g_map)
#> The optimized `gips` object.
#> 
#> Permutation:
#>  (1,2,3,4,5)
#> 
#> Log_posteriori:
#>  -0.5708509
#> 
#> Times more likely than starting permutation:
#>  953.095
#> 
#> The p-value of Likelihood-Ratio test:
#>  0.6008
#> 
#> The number of observations:
#>  13
#> 
#> The mean in the `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There are 12 degrees of freedom left.
#> 
#> n0 = number of cycles in a permutation + 1 (as the mean was estimated):
#>  2
#> 
#> The MLE estimator based on the found permutation does exist,
#> since the number of observations (13) is bigger than n0 (2),
#> and the projected covariance matrix is numerically positive definite.
#> 
#> The number of free parameters in the covariance matrix:
#>  3
#> 
#> BIC:
#>  104.0469
#> 
#> AIC:
#>  102.352
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  brute_force
#> 
#> The number of log_posteriori calls:
#>  67
#> 
#> Optimization time:
#>  0.04133129 secs
plot(g_map, type = "both", logarithmic_x = TRUE)
```
