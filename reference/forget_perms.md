# Forget the permutations for `gips` object optimized with `save_all_perms = TRUE`

Slim the `gips` object by forgetting the visited permutations from
`find_MAP(save_all_perms = TRUE)`.

## Usage

``` r
forget_perms(g)
```

## Arguments

- g:

  An object of class `gips`. A result of a
  `find_MAP(save_all_perms = TRUE)`.

## Value

The same object `g` as given, but without the visited permutation list.

## Details

For example, `perm_size = 150` and `max_iter = 150000` we checked
`forget_perms()` saves ~350 MB of RAM.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The `forget_perms()` is called on the output of
  `find_MAP(save_all_perms = TRUE)`.

## Examples

``` r
A <- matrix(rnorm(10 * 10), nrow = 10)
S <- t(A) %*% A
g <- gips(S, 13, was_mean_estimated = FALSE)
g_map <- find_MAP(g,
  max_iter = 10, optimizer = "Metropolis_Hastings",
  show_progress_bar = FALSE, save_all_perms = TRUE
)

object.size(g_map) # ~18 KB
#> 20672 bytes
g_map_slim <- forget_perms(g_map)
object.size(g_map_slim) # ~8 KB
#> 10680 bytes
```
