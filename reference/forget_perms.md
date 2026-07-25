# Forget the permutations for `gips` object optimized with `save_all_perms = TRUE`

Slim the `gips` object by forgetting the visited permutations from
`find_MAP(save_all_perms = TRUE)`.

## Usage

``` r
forget_perms(g)
```

## Arguments

- g:

  An object of class "gips"; a result of a
  `find_MAP(save_all_perms = TRUE)`.

## Value

Returns the same object `g` as given, but without the visited
permutation list.

## Details

For `perm_size = 150` and `max_iter = 150000` we checked it saves ~350
MB.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The `forget_perms()` is called on the output of
  `find_MAP(save_all_perms = TRUE)`.

## Examples

``` r
example_matrix <- matrix(rnorm(10 * 10), nrow = 10)
S <- t(example_matrix) %*% example_matrix
g <- gips(S, 13, was_mean_estimated = FALSE)
g_map <- find_MAP(g,
  max_iter = 10, optimizer = "Metropolis_Hastings",
  show_progress_bar = FALSE, save_all_perms = TRUE
)

object.size(g_map) # ~18 KB
#> 17688 bytes
g_map_slim <- forget_perms(g_map)
object.size(g_map_slim) # ~8 KB
#> 9008 bytes
```
