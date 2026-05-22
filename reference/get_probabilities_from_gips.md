# Extract probabilities for `gips` object optimized with `return_probabilities = TRUE`

After the `gips` object was optimized with the
`find_MAP(return_probabilities = TRUE)` function, then those calculated
probabilities can be extracted with this function.

## Usage

``` r
get_probabilities_from_gips(g)
```

## Arguments

- g:

  An object of class `gips`. A result of a
  `find_MAP(return_probabilities = TRUE)`.

## Value

Returns a numeric vector, calculated values of probabilities. Names
contain permutations this probabilities represent. For `gips` object
optimized with `find_MAP(return_probabilities = FALSE)`, it returns a
`NULL` object. It is sorted according to the probability.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The `get_probabilities_from_gips()` is called on the output of
  `find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`.

- [`vignette("Optimizers", package = "gips")`](https://przechoj.github.io/gips/articles/Optimizers.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Optimizers.html)) - A
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
