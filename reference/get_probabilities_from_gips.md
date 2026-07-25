# Extract probabilities for `gips` object optimized with `return_probabilities = TRUE`

After the `gips` object was optimized with
[`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
function with `return_probabilities = TRUE`, then those calculated
probabilities can be extracted with this function.

## Usage

``` r
get_probabilities_from_gips(g)
```

## Arguments

- g:

  An object of class "gips"; a result of a
  `find_MAP(return_probabilities = TRUE)`.

## Value

Returns a numeric vector, calculated values of probabilities. Names
contains permutations this probability represent. For `gips` object
optimized with `find_MAP(return_probabilities = FALSE)`, returns a
`NULL` object.

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
#>        ()     (1,2) 
#> 0.1649846 0.8350154 
```
