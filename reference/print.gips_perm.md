# Printing `gips_perm` object

Printing function for a `gips_perm` class.

## Usage

``` r
# S3 method for class 'gips_perm'
print(x, ...)
```

## Arguments

- x:

  An object of a `gips_perm` class.

- ...:

  Further arguments passed to
  [`permutations::print.cycle()`](https://robinhankin.github.io/permutations/reference/print.html).

## Value

Returns its argument invisibly, after printing it.

## Examples

``` r
g_perm <- gips_perm(permutations::as.cycle("(5,4)"), 5)
print(g_perm)
#> [1] (45)
```
