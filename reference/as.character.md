# Transform `gips_perm` object to character vector

Implementation of S3 method.

## Usage

``` r
# S3 method for class 'gips_perm'
as.character(x, ...)
```

## Arguments

- x:

  An object of a `gips_perm` class.

- ...:

  Further arguments passed to
  [`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html).

## Value

Returns an object of a `character` type.

## Methods (by class)

- `as.character(gips_perm)`:

## See also

[`permutations::as.character.cycle()`](https://robinhankin.github.io/permutations/reference/print.html)

## Examples

``` r
g_perm <- gips_perm(permutations::as.cycle("(5,4)"), 5)
as.character(g_perm)
#> [1] "(4,5)"
```
