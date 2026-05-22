# Printing `gips` object

Printing function for a `gips` class.

## Usage

``` r
# S3 method for class 'gips'
print(
  x,
  digits = 3,
  compare_to_original = TRUE,
  log_value = FALSE,
  oneline = FALSE,
  ...
)
```

## Arguments

- x:

  An object of a `gips` class.

- digits:

  The number of digits after the comma for a posteriori to be presented.
  It can be negative. By default, `Inf`. It is passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html).

- compare_to_original:

  A logical. Whether to print how many times more likely is the current
  permutation compared to:

  - the identity permutation `()` (for unoptimized `gips` object);

  - the starting permutation (for optimized `gips` object).

- log_value:

  A logical. Whether to print the logarithmic value. Default to `FALSE`.

- oneline:

  A logical. Whether to print in one or multiple lines. Default to
  `FALSE`.

- ...:

  The additional arguments passed to
  [`base::cat()`](https://rdrr.io/r/base/cat.html).

## Value

Returns an invisible `NULL`.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  The function that makes an optimized `gips` object out of the
  unoptimized one.

- [`compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md) -
  The function that prints the compared posteriories between any two
  permutations, not only compared to the starting one or id.

## Examples

``` r
S <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
g <- gips(S, 10, perm = "(12)")
print(g, digits = 4, oneline = TRUE)
#> The permutation (1,2): is 2.2801 times more likely than the () permutation.
```
