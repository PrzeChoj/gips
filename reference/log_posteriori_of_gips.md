# A log of a posteriori that the covariance matrix is invariant under permutation

More precisely, it is the logarithm of an unnormalized posterior
probability. It is the goal function for optimization algorithms in
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
function. The `perm_proposal` that maximizes this function is the
Maximum A Posteriori (MAP) Estimator.

## Usage

``` r
log_posteriori_of_gips(g)
```

## Arguments

- g:

  An object of a `gips_perm` class.

## Value

Returns a value of the logarithm of an unnormalized A Posteriori.

## Details

It is calculated using [formulas (33) and (27) from
references](https://arxiv.org/abs/2004.03503).

If `Inf` or `NaN` is reached, it produces a warning.

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- [`calculate_gamma_function()`](https://przechoj.github.io/gips/reference/calculate_gamma_function.md) -
  The function that calculates the value needed for
  `log_posteriori_of_gips()`.

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The functions that tries to optimize the `log_posteriori_of_gips`
  function.

- [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html) - A place
  to learn more about the math behind the `gips` package.

## Examples

``` r
# In the space with p = 2, there is only 2 permutations:
perm1 <- permutations::as.cycle(permutations::as.word(c(1, 2))) # (1)(2)
perm2 <- permutations::as.cycle(permutations::as.word(c(2, 1))) # (1,2)
S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
g1 <- gips(S1, 100, perm = perm1)
g2 <- gips(S1, 100, perm = perm2)
log_posteriori_of_gips(g1) # -136.6, this is the MAP Estimator
#> [1] -135.3433
log_posteriori_of_gips(g2) # -140.4
#> [1] -139.0036

exp(log_posteriori_of_gips(g1) - log_posteriori_of_gips(g2)) # 41.3
#> [1] 38.87122
# g1 is over 40 times more likely than g2.
# This is the expected outcome because S[1,1] significantly differs from S[2,2].

# ========================================================================

S2 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
g1 <- gips(S2, 100, perm = perm1)
g2 <- gips(S2, 100, perm = perm2)
log_posteriori_of_gips(g1) # -99.5
#> [1] -98.54173
log_posteriori_of_gips(g2) # -96.9, this is the MAP Estimator
#> [1] -96.00429

exp(log_posteriori_of_gips(g2) - log_posteriori_of_gips(g1)) # 12.7
#> [1] 12.64729
# g2 is over 12 times more likely than g1.
# This is the expected outcome because S[1,1] is very close to S[2,2].
```
