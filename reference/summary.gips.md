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

  3.  `start_permutation_log_posteriori` - the log of the A Posteriori
      value the start permutation has.

  4.  `times_more_likely_than_id` - how many times more likely the
      `start_permutation` is over the identity permutation, `()`. It can
      be less than 1, meaning the identity permutation is more likely.
      Remember that this number can become very large and overflow to
      `Inf` or small and underflow to 0.

  5.  `log_times_more_likely_than_id` - log of
      `times_more_likely_than_id`.

  6.  `likelihood_ratio_test_statistics`,
      `likelihood_ratio_test_p_value` - statistics and p-value of
      Likelihood Ratio test, where the H_0 is that the data was drawn
      from the normal distribution with Covariance matrix invariant
      under the given permutation. The p-value is calculated from the
      asymptotic distribution. Note that this is sensibly defined only
      for \\n \ge p\\.

  7.  `n0` - the minimum number of observations needed for the
      covariance matrix's maximum likelihood estimator (corresponding to
      a MAP) to exist. See **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html).

  8.  `S_matrix` - the underlying matrix. This matrix will be used in
      calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  9.  `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated. This value will be
      used in calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  10. `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was the output of
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function;

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`.

  11. `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
      See the **Hyperparameters** section of
      [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
      documentation.

  12. `n_parameters` - number of free parameters in the covariance
      matrix.

  13. `AIC`, `BIC` - output of
      [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      and
      [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      functions.

- For optimized `gips` object:

  1.  `optimized` - `TRUE`.

  2.  `found_permutation` - the permutation this `gips` represents. The
      visited permutation with the biggest A Posteriori value.

  3.  `found_permutation_log_posteriori` - the log of the A Posteriori
      value the found permutation has.

  4.  `start_permutation` - the original permutation this `gips`
      represented before optimization. It is the first visited
      permutation.

  5.  `start_permutation_log_posteriori` - the log of the A Posteriori
      value the start permutation has.

  6.  `times_more_likely_than_start` - how many more likely the
      `found_permutation` is over the `start_permutation`. It cannot be
      a number less than 1. Remember that this number can big and
      overflow to `Inf`.

  7.  `log_times_more_likely_than_start` - log of
      `times_more_likely_than_start`.

  8.  `likelihood_ratio_test_statistics`,
      `likelihood_ratio_test_p_value` - statistics and p-value of
      Likelihood Ratio test, where the H_0 is that the data was drawn
      from the normal distribution with Covariance matrix invariant
      under `found_permutation`. The p-value is calculated from the
      asymptotic distribution. Note that this is sensibly defined only
      for \\n \ge p\\.

  9.  `n0` - the minimal number of observations needed for the existence
      of the maximum likelihood estimator (corresponding to a MAP) of
      the covariance matrix (see **\\C\sigma\\ and `n0`** section in
      [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
      or in its [pkgdown
      page](https://przechoj.github.io/gips/articles/Theory.html)).

  10. `S_matrix` - the underlying matrix. This matrix will be used in
      calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  11. `number_of_observations` - the number of observations that were
      observed for the `S_matrix` to be calculated. This value will be
      used in calculations of the posteriori value in
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

  12. `was_mean_estimated` - given by the user while creating the `gips`
      object:

      - `TRUE` means the `S` parameter was output of the
        [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function;

      - `FALSE` means the `S` parameter was calculated with
        `S = t(X) %*% X / number_of_observations`.

  13. `delta`, `D_matrix` - the hyperparameters of the Bayesian method.
      See the **Hyperparameters** section of
      [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
      documentation.

  14. `n_parameters` - number of free parameters in the covariance
      matrix.

  15. `AIC`, `BIC` - output of
      [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      and
      [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
      functions.

  16. `optimization_algorithm_used` - all used optimization algorithms
      in order (one could start optimization with "MH", and then do an
      "HC").

  17. `did_converge` - a boolean, did the last used algorithm converge.

  18. `number_of_log_posteriori_calls` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called during the optimization.

  19. `whole_optimization_time` - how long was the optimization process;
      the sum of all optimization times (when there were multiple).

  20. `log_posteriori_calls_after_best` - how many times was the
      [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
      function called after the `found_permutation`; in other words, how
      long ago could the optimization be stopped and have the same
      result. If this value is small, consider running
      [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
      again with `optimizer = "continue"`. For `optimizer = "BF"`, it is
      `NULL`.

  21. `acceptance_rate` - only interesting for `optimizer = "MH"`. How
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
    1.1, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.1, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.1, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.1, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.1, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.1
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)
unclass(summary(g))
#> $optimized
#> [1] FALSE
#> 
#> $start_permutation
#> [1] ()
#> 
#> $start_permutation_log_posteriori
#> [1] -31.54222
#> 
#> $times_more_likely_than_id
#> [1] 1
#> 
#> $log_times_more_likely_than_id
#> [1] 0
#> 
#> $likelihood_ratio_test_statistics
#> [1] 0
#> 
#> $likelihood_ratio_test_p_value
#> NULL
#> 
#> $n0
#> [1] 7
#> 
#> $S_matrix
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 1.0494842 1.1827957 1.0811047 0.6491538 0.5910505 0.7578313
#> [2,] 1.1827957 1.5722816 1.1559233 0.5854755 0.2377610 0.7524728
#> [3,] 1.0811047 1.1559233 1.5841572 0.9470600 0.7940677 0.5350399
#> [4,] 0.6491538 0.5854755 0.9470600 1.0954129 0.9148409 0.6746779
#> [5,] 0.5910505 0.2377610 0.7940677 0.9148409 1.5101372 0.7777238
#> [6,] 0.7578313 0.7524728 0.5350399 0.6746779 0.7777238 1.0100209
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
#> [1,] 1.303582 0.000000 0.000000 0.000000 0.000000 0.000000
#> [2,] 0.000000 1.303582 0.000000 0.000000 0.000000 0.000000
#> [3,] 0.000000 0.000000 1.303582 0.000000 0.000000 0.000000
#> [4,] 0.000000 0.000000 0.000000 1.303582 0.000000 0.000000
#> [5,] 0.000000 0.000000 0.000000 0.000000 1.303582 0.000000
#> [6,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.303582
#> 
#> $n_parameters
#> [1] 21
#> 
#> $AIC
#> [1] 163.442
#> 
#> $BIC
#> [1] 175.3059
#> 

g_map <- find_MAP(g, max_iter = 10, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
unclass(summary(g_map))
#> $optimized
#> [1] TRUE
#> 
#> $found_permutation
#> [1] (12)(45)
#> 
#> $found_permutation_log_posteriori
#> [1] -28.4427
#> 
#> $start_permutation
#> [1] ()
#> 
#> $start_permutation_log_posteriori
#> [1] -31.54222
#> 
#> $times_more_likely_than_start
#> [1] 22.18745
#> 
#> $log_times_more_likely_than_start
#> [1] 3.099527
#> 
#> $likelihood_ratio_test_statistics
#> [1] 17.4246
#> 
#> $likelihood_ratio_test_p_value
#> [1] 0.0259792
#> 
#> $n0
#> [1] 5
#> 
#> $S_matrix
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 1.0494842 1.1827957 1.0811047 0.6491538 0.5910505 0.7578313
#> [2,] 1.1827957 1.5722816 1.1559233 0.5854755 0.2377610 0.7524728
#> [3,] 1.0811047 1.1559233 1.5841572 0.9470600 0.7940677 0.5350399
#> [4,] 0.6491538 0.5854755 0.9470600 1.0954129 0.9148409 0.6746779
#> [5,] 0.5910505 0.2377610 0.7940677 0.9148409 1.5101372 0.7777238
#> [6,] 0.7578313 0.7524728 0.5350399 0.6746779 0.7777238 1.0100209
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
#> [1,] 1.303582 0.000000 0.000000 0.000000 0.000000 0.000000
#> [2,] 0.000000 1.303582 0.000000 0.000000 0.000000 0.000000
#> [3,] 0.000000 0.000000 1.303582 0.000000 0.000000 0.000000
#> [4,] 0.000000 0.000000 0.000000 1.303582 0.000000 0.000000
#> [5,] 0.000000 0.000000 0.000000 0.000000 1.303582 0.000000
#> [6,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.303582
#> 
#> $n_parameters
#> [1] 13
#> 
#> $AIC
#> [1] 164.8666
#> 
#> $BIC
#> [1] 172.211
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
#> Time difference of 0.01308608 secs
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
#>  (3,6)(4,5)
#> 
#> Log_posteriori:
#>  -25.37572
#> 
#> Times more likely than starting permutation:
#>  476.515
#> 
#> The p-value of Likelihood-Ratio test:
#>  0.2636
#> 
#> The number of observations:
#>  13
#> 
#> The mean in the `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There are 12 degrees of freedom left.
#> 
#> n0:
#>  5
#> 
#> The number of observations is bigger than n0 for this permutation,
#> so the gips model based on the found permutation does exist.
#> 
#> The number of free parameters in the covariance matrix:
#>  13
#> 
#> BIC:
#>  164.8062
#> 
#> AIC:
#>  157.4618
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  hill_climbing did converge
#> 
#> The number of log_posteriori calls:
#>  46
#> 
#> Optimization time:
#>  0.06210971 secs
#> 
#> Log_posteriori calls after the found permutation:
#>  17
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
#> The current permutation is id, so Likelihood-Ratio test cannot be performed (there is nothing to compare)
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
