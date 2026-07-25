# Extract the Log-Likelihood for `gips` class

Calculates Log-Likelihood of the sample based on the `gips` object.

## Usage

``` r
# S3 method for class 'gips'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `gips`. Usually, a result of a
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- ...:

  Further arguments will be ignored.

## Value

Log-Likelihood of the sample. Object of class `logLik`.

Possible failure situations:

- When the multivariate normal model does not exist
  (`number_of_observations < n0`), it returns `NULL`.

- When the multivariate normal model cannot be reasonably approximated
  (output of
  [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
  is singular), it returns `-Inf`.

In both failure situations, it shows a warning. More information can be
found in the **Existence of likelihood** section below.

## Details

This will always be the biggest for `perm = "()"` (provided that
`p <= n`).

If the found permutation still requires more parameters than `n`, the
likelihood does not exist; thus the function returns `NULL`.

If the `projected_cov` (output of
[`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md))
is close to singular, the `NA` is returned.

## Existence of likelihood

We only consider the non-degenerate multivariate normal model. In the
`gips` context, such a model exists only when the number of observations
is bigger or equal to `n0`. To get `n0` for the `gips` object `g`, call
`summary(g)$n0`.

See examples where the `g_n_too_small` had too small
`number_of_observations` to have likelihood. After the optimization, the
likelihood did exist.

For more information, refer to **\\C\_\sigma\\ and `n0`** section in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html).

## Calculation details

For more details and used formulas, see the **Information Criterion -
AIC and BIC** section in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html).

## See also

- [`logLik()`](https://rdrr.io/r/stats/logLik.html) - Generic function
  this `logLik.gips()` extends.

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  Usually, the `logLik.gips()` is called on the output of
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md),
  [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md) -
  Often, one is more interested in an Information Criterion AIC or BIC.

- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md) -
  One can get `n0` by calling `summary(g)$n0`. To see why one may be
  interested in `n0`, see the **Existence of likelihood** section above.

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  Project the known matrix onto the found permutations space. It is
  mentioned in the **Calculation details** section above.

## Examples

``` r
S <- matrix(c(
  5.15, 2.05, 3.60, 1.99,
  2.05, 5.09, 2.03, 3.57,
  3.60, 2.03, 5.21, 1.97,
  1.99, 3.57, 1.97, 5.13
), nrow = 4)
g <- gips(S, 5)
logLik(g) # -32.67048
#> 'log Lik.' -32.67048 (df=10)
# For perm = "()", which is default, there is p + choose(p, 2) degrees of freedom

g_map <- find_MAP(g, optimizer = "brute_force")
#> ================================================================================
logLik(g_map) # -32.6722 # this will always be smaller than `logLik(gips(S, n, perm = ""))`
#> 'log Lik.' -32.6722 (df=3)

g_n_too_small <- gips(S, number_of_observations = 4)
logLik(g_n_too_small) # NULL # the likelihood does not exists
#> Warning: The likelihood is not defined for this `gips`.
#> ✖ The n = 4 is smaller than the minimum required n0 = 5. For more information, see section **Existence of likelihood** in documentation `?logLik.gips` or its [pkgdown page](https://przechoj.github.io/gips/reference/logLik.gips.html).
#> NULL
summary(g_n_too_small)$n0 # 5, but we set number_of_observations = 4, which is smaller
#> [1] 5

g_MAP <- find_MAP(g_n_too_small)
#> You used the default value of the 'optimizer' argument in `find_MAP()`.
#> ℹ The 'optimizer = NA' was automatically changed to 'optimizer = "BF"'.
#> ================================================================================
logLik(g_MAP) # -24.94048, this is no longer NULL
#> 'log Lik.' -24.94048 (df=3)
summary(g_MAP)$n0 # 2
#> [1] 2
```
