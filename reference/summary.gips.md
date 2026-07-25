# Summarizing the gips object

`summary` method for `gips` class.

## Usage

``` r
# S3 method for class 'gips'
summary(object, ...)

# S3 method for class 'summary.gips'
print(x, ...)
```

## Arguments

- object:

  An object of class `gips`. Usually, a result of a
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- ...:

  Further arguments passed to or from other methods.

- x:

  An object of class `summary.gips` to be printed

## Value

The function `summary.gips()` computes and returns a list of summary
statistics of the given `gips` object. Those are:

- For unoptimized `gips` object:

  1.  `optimized` - `FALSE`.

  2.  `start_permutation` - the permutation this `gips` represents.

  3.  `start_permutation_log_posteriori` - the log of the a posteriori
      value the start permutation has.

  4.  `times_more_likely_than_id` - how many more likely the
      `start_permutation` is over the identity permutation, `()`. It can
      be less than 1, meaning the identity permutation is more likely.
      Remember that this number can big and overflow to `Inf` or small
      and underflow to 0.

  5.  `n0` - the minimum number of observations needed for the
      covariance matrix's maximum likelihood estimator (corresponding to
      a MAP) to exist. See **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html).

  6.  `S_matrix` - the underlying matrix. This matrix will be used in
      calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  7.  `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated. This value will be
      used in calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  8.  `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was the output of
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function;

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`.

  9.  `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
      See the **Hyperparameters** section of
      [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
      documentation.

  10. `AIC`, `BIC` - output of
      [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      and
      [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      functions.

  11. `n_parameters` - number of free parameters in the covariance
      matrix.

- For optimized `gips` object:

  1.  `optimized` - `TRUE`.

  2.  `found_permutation` - the permutation this `gips` represents. The
      visited permutation with the biggest a posteriori value.

  3.  `found_permutation_log_posteriori` - the log of the a posteriori
      value the found permutation has.

  4.  `start_permutation` - the original permutation this `gips`
      represented before optimization. It is the first visited
      permutation.

  5.  `start_permutation_log_posteriori` - the log of the a posteriori
      value the start permutation has.

  6.  `times_more_likely_than_start` - how many more likely the
      `found_permutation` is over the `start_permutation`. It cannot be
      a number less than 1. Remember that this number can big and
      overflow to `Inf`.

  7.  `n0` - the minimal number of observations needed for the existence
      of the maximum likelihood estimator (corresponding to a MAP) of
      the covariance matrix (see **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html)).

  8.  `S_matrix` - the underlying matrix. This matrix will be used in
      calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  9.  `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated. This value will be
      used in calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  10. `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was output of the
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function;

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`.

  11. `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
      See the **Hyperparameters** section of
      [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
      documentation.

  12. `AIC`, `BIC` - output of
      [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      and
      [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      functions.

  13. `n_parameters` - number of free parameters in the covariance
      matrix.

  14. `optimization_algorithm_used` - all used optimization algorithms
      in order (one could start optimization with "MH", and then do an
      "HC").

  15. `did_converge` - a boolean, did the last used algorithm converge.

  16. `number_of_log_posteriori_calls` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called during the optimization.

  17. `whole_optimization_time` - how long was the optimization process;
      the sum of all optimization times (when there were multiple).

  18. `log_posteriori_calls_after_best` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called after the `found_permutation`; in other words, how
      long ago could the optimization be stopped and have the same
      result. If this value is small, consider running
      [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
      again with `optimizer = "continue"`. For `optimizer = "BF"`, it is
      `NULL`.

  19. `acceptance_rate` - only interesting for `optimizer = "MH"`. How
      often was the algorithm accepting the change of permutation in an
      iteration.

The function `print.summary.gips()` returns an invisible `NULL`.

## Methods (by generic)

- `print(summary.gips)`: Printing method for class `summary.gips`.
  Prints every interesting information in a form pleasant for humans.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  Usually, the `summary.gips()` is called on the output of
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  Calculate the likelihood of a permutation.

- [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md),
  [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md) -
  Calculate Akaike's or Bayesian Information Criterion

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  Project the known matrix on the found permutations space.

## Examples

``` r
require("MASS") # for mvrnorm()

perm_size <- 6
mu <- runif(6, -10, 10) # Assume we don't know the mean
sigma_matrix <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)

g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
unclass(summary(g_map))
#> $optimized
#> [1] TRUE
#> 
#> $found_permutation
#> [1] ()
#> 
#> $found_permutation_log_posteriori
#> [1] -4.536729
#> 
#> $start_permutation
#> [1] ()
#> 
#> $start_permutation_log_posteriori
#> [1] -4.536729
#> 
#> $times_more_likely_than_start
#> [1] 1
#> 
#> $log_times_more_likely_than_start
#> [1] 0
#> 
#> $n0
#> [1] 7
#> 
#> $S_matrix
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 0.5231376 0.3721420 0.4435844 0.7453107 0.8963063 0.8248638
#> [2,] 0.3721420 0.5271158 0.6757445 0.9053160 0.7503422 0.6017135
#> [3,] 0.4435844 0.6757445 1.0217525 1.2731260 1.0409659 0.6949579
#> [4,] 0.7453107 0.9053160 1.2731260 1.9798133 1.8198079 1.4519979
#> [5,] 0.8963063 0.7503422 1.0409659 1.8198079 1.9657720 1.6751482
#> [6,] 0.8248638 0.6017135 0.6949579 1.4519979 1.6751482 1.5819038
#> 
#> $number_of_observations
#> [1] 13
#> 
#> $was_mean_estimated
#> [1] TRUE
#> 
#> $delta
#> [1] 3
#> 
#> $D_matrix
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#> [1,] 1.266582 0.000000 0.000000 0.000000 0.000000 0.000000
#> [2,] 0.000000 1.266582 0.000000 0.000000 0.000000 0.000000
#> [3,] 0.000000 0.000000 1.266582 0.000000 0.000000 0.000000
#> [4,] 0.000000 0.000000 0.000000 1.266582 0.000000 0.000000
#> [5,] 0.000000 0.000000 0.000000 0.000000 1.266582 0.000000
#> [6,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.266582
#> 
#> $n_parameters
#> [1] 21
#> 
#> $AIC
#> [1] -653.0269
#> 
#> $BIC
#> [1] -641.163
#> 
#> $optimization_algorithm_used
#> [1] "Metropolis_Hastings"
#> 
#> $did_converge
#> NULL
#> 
#> $number_of_log_posteriori_calls
#> [1] 10
#> 
#> $whole_optimization_time
#> Time difference of 0.01458597 secs
#> 
#> $log_posteriori_calls_after_best
#> [1] 9
#> 
#> $acceptance_rate
#> [1] 0.3
#> 

g_map2 <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "hill_climbing")
summary(g_map2)
#> The optimized `gips` object.
#> 
#> Permutation:
#>  ()
#> 
#> Log_posteriori:
#>  -4.536729
#> 
#> Times more likely than starting permutation:
#>  1
#> 
#> The number of observations:
#>  13
#> 
#> The mean in the `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There are 12 degrees of freedom left.
#> 
#> n0:
#>  7
#> 
#> The number of observations is bigger than n0 for this permutation,
#> so the gips model based on the found permutation does exist.
#> 
#> The number of free parameters in the covariance matrix:
#>  21
#> 
#> BIC:
#>  -641.163
#> 
#> AIC:
#>  -653.0269
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  hill_climbing did converge
#> 
#> The number of log_posteriori calls:
#>  16
#> 
#> Optimization time:
#>  0.0247817 secs
#> 
#> Log_posteriori calls after the found permutation:
#>  15
# ================================================================================
S <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
g <- gips(S, 10)
print(summary(g))
#> The unoptimized `gips` object.
#> 
#> Permutation:
#>  ()
#> 
#> Log_posteriori:
#>  -15.4837
#> 
#> The number of observations:
#>  10
#> 
#> The mean in the `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There are 9 degrees of freedom left.
#> 
#> n0:
#>  3
#> 
#> The number of observations is bigger than n0 for this permutation,
#> so the gips model based on the found permutation does exist.
#> 
#> The number of free parameters in the covariance matrix:
#>  3
#> 
#> BIC:
#>  63.02608
#> 
#> AIC:
#>  62.11833
```
