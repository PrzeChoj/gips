# Summarizing the gips object

`summary` method for class "gips".

## Usage

``` r
# S3 method for class 'gips'
summary(object, ...)

# S3 method for class 'summary.gips'
print(x, ...)
```

## Arguments

- object:

  An object of class "gips"; is usually a result of a
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- ...:

  Further arguments passed to or from other methods.

- x:

  An object of class "summary.gips" to be printed

## Value

The function `summary.gips` computes and returns a list of summary
statistics of the given `gips` object. Those are:

- For unoptimized `gips` object:

  1.  `optimized` - `FALSE`

  2.  `start_permutation` - the permutation this `gips` represents

  3.  `start_permutation_log_posteriori` - the log of the a posteriori
      value the start permutation has

  4.  `times_more_likely_than_id` - how many more likely the
      `start_permutation` is over the identity permutation, `()`. It can
      be a number less than 1, which means the identity permutation,
      `()`, is more likely. Keep in mind this number can be really big
      and can be overflowed to `Inf`

  5.  `n0` - the minimal number of observations needed for existence of
      the maximum likelihood estimator (corresponding to a MAP) of the
      covariance matrix (see **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html)).

  6.  `S_matrix` - the underlying matrix; this is used to calculate the
      posteriori value

  7.  `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated; this is used to
      calculate the posteriori value

  8.  `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was output of
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`

  9.  `delta`, `D_matrix` - the parameters of the Bayesian method

- For optimized `gips` object:

  1.  `optimized` - `TRUE`

  2.  `found_permutation` - the permutation this `gips` represents; the
      visited permutation with the biggest a posteriori value

  3.  `found_permutation_log_posteriori` - the log of the a posteriori
      value the found permutation have

  4.  `start_permutation` - the original permutation this `gips`
      represented before optimization; the first visited permutation

  5.  `start_permutation_log_posteriori` - the log of the a posteriori
      value the start permutation has

  6.  `times_more_likely_than_start` - how many more likely the
      `found_permutation` is over the `start_permutation`. It cannot be
      a number less than 1. Keep in mind this number can be really big
      and can be overflowed to `Inf`

  7.  `n0` - the minimal number of observations needed for existence of
      the maximum likelihood estimator (corresponding to a MAP) of the
      covariance matrix (see **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html)).

  8.  `S_matrix` - the underlying matrix; this is used to calculate the
      posteriori value

  9.  `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated; this is used to
      calculate the posteriori value

  10. `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was output of
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`

  11. `delta`, `D_matrix` - the parameters of the Bayesian method

  12. `optimization_algorithm_used` - all used optimization algorithms
      in order (one could start optimization with "MH", and then do an
      "HC")

  13. `did_converge` - a boolean, did the last used algorithm converge

  14. `number_of_log_posteriori_calls` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called during the optimization

  15. `whole_optimization_time` - how long was the optimization process;
      the sum of all optimization times (when there were multiple)

  16. `log_posteriori_calls_after_best` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called after the `found_permutation`; in other words, how
      long ago could the optimization be stopped and have the same
      result; if this value is small, consider running
      [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
      one more time with `optimizer = "continue"`. For
      `optimizer = "BF"`, it is `NULL`

  17. `acceptance_rate` - only interesting for `optimizer = "MH"`; how
      often was the algorithm accepting the change of permutation in an
      iteration

`print.summary.gips` returns an invisible `NULL`.

## Methods (by generic)

- `print(summary.gips)`: Printing method for class "summary.gips".
  Prints every interesting information in a pleasant, human readable
  form

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  Usually, the `summary.gips()` is called on the output of
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  The function that calculates the likelihood of a permutation.

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  The function that can project the known matrix of the found
  permutations space.

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
#> [1] (25)(36)
#> 
#> $found_permutation_log_posteriori
#> [1] -5.171758
#> 
#> $start_permutation
#> [1] ()
#> 
#> $start_permutation_log_posteriori
#> [1] -9.436118
#> 
#> $times_more_likely_than_start
#> [1] 71.11942
#> 
#> $n0
#> [1] 5
#> 
#> $S_matrix
#>            [,1]        [,2]      [,3]       [,4]        [,5]      [,6]
#> [1,] 0.68035238  0.43179370 0.2145354 0.08460627  0.33316494 0.5504233
#> [2,] 0.43179370  0.81291728 0.5596823 0.33203041 -0.04909318 0.2041418
#> [3,] 0.21453537  0.55968234 0.5546893 0.45125412  0.10610714 0.1111001
#> [4,] 0.08460627  0.33203041 0.4512541 0.70575290  0.45832876 0.3391051
#> [5,] 0.33316494 -0.04909318 0.1061071 0.45832876  0.84058688 0.6853866
#> [6,] 0.55042328  0.20414176 0.1111001 0.33910506  0.68538657 0.7784282
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
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    0
#> [2,]    0    1    0    0    0    0
#> [3,]    0    0    1    0    0    0
#> [4,]    0    0    0    1    0    0
#> [5,]    0    0    0    0    1    0
#> [6,]    0    0    0    0    0    1
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
#> Time difference of 0.01854467 secs
#> 
#> $log_posteriori_calls_after_best
#> [1] 0
#> 
#> $acceptance_rate
#> [1] 0.2
#> 

g_map2 <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "hill_climbing")
summary(g_map2)
#> The optimized `gips` object.
#> 
#> Permutation:
#>  (1,4)
#> 
#> Log_posteriori:
#>  -7.453653
#> 
#> Times more likely than starting permutation:
#> 7.26061831298498
#> 
#> Number of observations:
#>  13
#> 
#> The mean in `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There is 12 degrees of freedom left.
#> 
#> n0:
#>  6
#> 
#> Number of observations is bigger than n0 for this permutaion,
#> so the gips model based on the found permutation does exist.
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  hill_climbing did converge
#> 
#> Number of log_posteriori calls:
#>  31
#> 
#> Optimization time:
#>  0.04063201 secs
#> 
#> Log_posteriori calls after the found permutation:
#>  27
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
#>  -16.56396
#> 
#> Times more likely than identity permutation:
#> 1
#> 
#> Number of observations:
#>  10
#> 
#> The mean in `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There is 9 degrees of freedom left.
#> 
#> n0:
#>  3
#> 
#> Number of observations is bigger than n0 for this permutaion,
#> so the gips model based on the found permutation does exist.
```
