# The constructor of a `gips` class.

Create a `gips` object. This object will consist of data and all other
information needed to find the most likely invariant permutation. The
optimization itself will not be performed. One must call the
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
function to do it. See examples below.

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

  A matrix; estimated covariance matrix. When `Z` is the observed data:

  - if one does not know the theoretical mean and has to estimate it
    with the observed mean, use `S = cov(Z)`, and leave parameter
    `was_mean_estimated = TRUE` as default.

  - if one know the theoretical mean is 0, use
    `S = (t(Z) %*% Z) / number_of_observations`, and set parameter
    `was_mean_estimated = FALSE`;

- number_of_observations:

  A number of data points that `S` is based on.

- delta:

  A number, hyper-parameter of a Bayesian model. Has to be bigger
  than 2. See **Hyperparameters** section bellow.

- D_matrix:

  A symmetric, positive-definite matrix of the same size as `S`.
  Hyper-parameter of a Bayesian model. When `NULL`, the identity matrix
  is taken. See **Hyperparameters** section bellow.

- was_mean_estimated:

  A boolean.

  - Set `TRUE` (default) when your `S` parameter is a result of a
    [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function.

  - Set FALSE when your `S` parameter is a result of a
    `(t(Z) %*% Z) / number_of_observations` calculation.

- perm:

  An optional permutation to be the base for the `gips` object. Can be
  of a `gips_perm` or a `permutation` class, or anything the function
  [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html)
  can handle.

- list_of_gips_perm:

  A list with a single element of a `gips_perm` class. The base object
  for the `gips` object.

- optimization_info:

  For internal use only. `NULL` or the list with information about the
  optimization process.

- g:

  Object to be checked whether it is proper object of a `gips` class.

## Value

`gips()` returns an object of a `gips` class after the safety checks.

`new_gips()` returns an object of a `gips` class without the safety
checks.

`validate_gips()` returns its argument unchanged. If the argument is not
a proper element of a `gips` class, it produces an error.

## Functions

- `new_gips()`: Constructor. Only intended for low-level use.

- `validate_gips()`: Validator. Only intended for low-level use.

## Methods for a `gips` class

- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md)

- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md)

- [`print.gips()`](https://przechoj.github.io/gips/reference/print.gips.md)

## Hyperparameters

In the Bayesian model, the prior distribution for the covariance matrix
is a generalized case of [Wishart
distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

For brief introduction, see **Bayesian model selection** section in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or in its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html)).

## See also

- [`stats::cov()`](https://rdrr.io/r/stats/cor.html) - The `S` parameter
  is most of the time an estimated covariance matrix, so a result of the
  [`cov()`](https://rdrr.io/r/stats/cor.html) function. For more
  information, see [Wikipedia - Estimation of covariance
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
#> The permutation (1,4,2,5,3)
#>  - was found after 120 log_posteriori calculations
#>  - is 19478.229036143 times more likely than the starting, () permutation.

summary(g_map)
#> The optimized `gips` object.
#> 
#> Permutation:
#>  (1,4,2,5,3)
#> 
#> Log_posteriori:
#>  -14.2724
#> 
#> Times more likely than starting permutation:
#> 19478.229036143
#> 
#> Number of observations:
#>  13
#> 
#> The mean in `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There is 12 degrees of freedom left.
#> 
#> n0:
#>  2
#> 
#> Number of observations is bigger than n0 for this permutaion,
#> so the gips model based on the found permutation does exist.
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  brute_force
#> 
#> Number of log_posteriori calls:
#>  120
#> 
#> Optimization time:
#>  0.1820033 secs

if (require("graphics")) {
  plot(g_map, type = "both", logarithmic_x = TRUE)
}
```
