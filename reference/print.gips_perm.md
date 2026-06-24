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

  Further arguments (currently ignored).

## Value

An invisible `NULL`.

## Examples

``` r
gperm <- gips_perm("(5,4)", 5)
print(gperm)
#> [1] (45)
```
