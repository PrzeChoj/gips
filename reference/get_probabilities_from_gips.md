# Extract probabilities for `gips` object optimized with `return_probabilities = TRUE`

Read probabilities previously calculated by
`find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`. This
function does not calculate missing probabilities. BF gives the full
posterior over cyclic groups; MH gives the posterior conditional on the
distinct visited groups. See
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

## Usage

``` r
get_probabilities_from_gips(g)
```

## Arguments

- g:

  An object of class `gips`. A result of a
  `find_MAP(return_probabilities = TRUE)`.

## Value

A numeric vector of probabilities sorted in decreasing order, named by
representative generators of cyclic groups. Returns `NULL` with a
message if the optimized object has no saved probabilities.

## Details

Saved visit-frequency estimates from older MH objects are returned
unchanged with a warning; run a new optimization for weight-based
results.

## Multi-sample

For multi-sample `gips` objects, these are probabilities of the shared
permutation symmetry, calculated from the combined posterior used by
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The `get_probabilities_from_gips()` is called on the output of
  `find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`.

- [`vignette("Optimizers", package = "gips")`](https://przechoj.github.io/gips/articles/Optimizers.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Optimizers.html) - A
  place to learn more about the available optimizers.

## Examples

``` r
g <- gips(matrix(c(1, 0.5, 0.5, 1.3), nrow = 2), 13, was_mean_estimated = FALSE)
g_map <- find_MAP(g,
  optimizer = "BF", show_progress_bar = FALSE,
  return_probabilities = TRUE, save_all_perms = TRUE
)

get_probabilities_from_gips(g_map)
#>     (1,2)        () 
#> 0.8170484 0.1829516 
```
