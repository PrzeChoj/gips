# Transform the `gips_perm` object to a character vector

Implementation of the S3 method.

## Usage

``` r
# S3 method for class 'gips_perm'
as.character(x, ...)
```

## Arguments

- x:

  An object of a `gips_perm` class.

- ...:

  Further arguments (currently ignored).

## Value

Returns an object of a `character` type.

## See also

- [`as.character.gips()`](https://przechoj.github.io/gips/reference/as.character.gips.md) -
  The underlying `gips_perm` of the `gips` object is passed to
  `as.character.gips_perm()`.

- [`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html) -
  The underlying permutation of the `gips` object is passed to
  [`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html).

## Examples

``` r
g_perm <- gips_perm("(5,4)", 5)
as.character(g_perm)
#> [1] "(4,5)"
```
