# Akaike's An Information Criterion for `gips` class

Akaike's An Information Criterion for `gips` class

## Usage

``` r
# S3 method for class 'gips'
AIC(object, ..., k = 2)

# S3 method for class 'gips'
BIC(object, ...)
```

## Arguments

- object:

  An object of class `gips`. Usually, a result of a
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- ...:

  Further arguments will be ignored.

- k:

  Numeric, the *penalty* per parameter to be used. The default `k = 2`
  is the classical AIC.

## Value

`AIC.gips()` returns calculated Akaike's An Information Criterion

When the multivariate normal model does not exist
(`number_of_observations < n0`), it returns `NULL`. When the
multivariate normal model cannot be reasonably approximated (output of
[`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
is singular), it returns `Inf`.

In both failure situations, shows a warning. More information can be
found in the **Existence of likelihood** section of
[`logLik.gips()`](https://przechoj.github.io/gips/reference/logLik.gips.md).

`BIC.gips()` returns calculated Schwarz's Bayesian Information
Criterion.

## Functions

- `BIC(gips)`: Schwarz's Bayesian Information Criterion

## Calculation details

For more details and used formulas, see the **Information Criterion -
AIC and BIC** section in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html).

## See also

- [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) - Generic functions this
  `AIC.gips()` and `BIC.gips()` extend.

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  Usually, the `AIC.gips()` and `BIC.gips()` are called on the output of
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- [`logLik.gips()`](https://przechoj.github.io/gips/reference/logLik.gips.md) -
  Calculates the log-likelihood for the `gips` object. An important part
  of the Information Criteria.

## Examples

``` r
S <- matrix(c(
  5.15, 2.05, 3.10, 1.99,
  2.05, 5.09, 2.03, 3.07,
  3.10, 2.03, 5.21, 1.97,
  1.99, 3.07, 1.97, 5.13
), nrow = 4)
g <- gips(S, 14)
g_map <- find_MAP(g, optimizer = "brute_force")
#> ================================================================================

AIC(g) # 238
#> [1] 237.6098
AIC(g_map) # 224 < 238, so g_map is better than g according to AIC
#> [1] 223.6188
# ================================================================================
BIC(g) # 244
#> [1] 244.0004
BIC(g_map) # 226 < 244, so g_map is better than g according to BIC
#> [1] 225.536
```
