# Transform the `gips` object to a character vector

Implementation of the S3 method.

## Usage

``` r
# S3 method for class 'gips'
as.character(x, ...)
```

## Arguments

- x:

  An object of a `gips` class.

- ...:

  Further arguments (currently ignored).

## Value

Returns an object of a `character` type.

## See also

- [`as.character.gips_perm()`](https://przechoj.github.io/gips/reference/as.character.gips_perm.md) -
  The underlying `gips_perm` of the `gips` object is passed to
  [`as.character.gips_perm()`](https://przechoj.github.io/gips/reference/as.character.gips_perm.md).

- [`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html) -
  The underlying permutation of the `gips` object is passed to
  [`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html).

## Examples

``` r
A <- matrix(rnorm(4 * 4), nrow = 4)
S <- t(A) %*% A
g <- gips(S, 14, perm = "(123)")
as.character(g)
#> [1] "(1,2,3)"
```
