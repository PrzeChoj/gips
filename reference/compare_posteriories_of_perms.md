# Compare the posteriori probabilities of 2 permutations

Check which permutation is more likely and how much more likely.

## Usage

``` r
compare_posteriories_of_perms(
  perm1,
  perm2 = "()",
  S = NULL,
  number_of_observations = NULL,
  delta = 3,
  D_matrix = NULL,
  was_mean_estimated = TRUE,
  print_output = TRUE,
  digits = 3
)

compare_log_posteriories_of_perms(
  perm1,
  perm2 = "()",
  S = NULL,
  number_of_observations = NULL,
  delta = 3,
  D_matrix = NULL,
  was_mean_estimated = TRUE,
  print_output = TRUE,
  digits = 3
)
```

## Arguments

- perm1, perm2:

  Permutations to compare. How many times `perm1` is more likely than
  `perm2`? Those can be provided as the `gips` objects, the `gips_perm`
  objects, or anything that can be used as the `x` parameter in the
  [`gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md)
  function. They do not have to be of the same class.

- S, number_of_observations, delta, D_matrix, was_mean_estimated:

  The same parameters as in the
  [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
  function. If at least one of `perm1` or `perm2` is a `gips` object,
  they are overwritten with those from the `gips` object.

- print_output:

  A boolean. When `TRUE` (default), the computed value will be printed
  with additional text and returned invisibly. When `FALSE`, the
  computed value will be returned visibly.

- digits:

  Integer. Only used when `print_output = TRUE`. The number of digits
  after the comma to print. It can be negative, can be `+Inf`. It is
  passed to [`base::round()`](https://rdrr.io/r/base/Round.html).

## Value

The function `compare_posteriories_of_perms()` returns the value of how
many times the `perm1` is more likely than `perm2`.

The function `compare_log_posteriories_of_perms()` returns the logarithm
of how many times the `perm1` is more likely than `perm2`.

## Functions

- `compare_log_posteriories_of_perms()`: More stable, logarithmic
  version of `compare_posteriories_of_perms()`. The natural logarithm is
  used.

## See also

- [`print.gips()`](https://przechoj.github.io/gips/reference/print.gips.md) -
  The function that prints the posterior of the optimized `gips` object
  compared to the starting permutation.

- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md) -
  The function that calculates the posterior of the optimized `gips`
  object compared to the starting permutation.

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The function that finds the permutation that maximizes
  [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  The function this `compare_posteriories_of_perms()` calls underneath.

## Examples

``` r
require("MASS") # for mvrnorm()
#> Loading required package: MASS

perm_size <- 6
mu <- runif(6, -10, 10) # Assume we don't know the mean
sigma_matrix <- matrix(
  data = c(
    1.05, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.05, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.05, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.05, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.05, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.05
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)
g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")

compare_posteriories_of_perms(g_map, g, print_output = TRUE)
#> The permutation () is 1 times more likely than the () permutation.
exp(compare_log_posteriories_of_perms(g_map, g, print_output = FALSE))
#> [1] 1
```
