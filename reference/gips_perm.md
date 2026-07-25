# Permutation object

Create permutation objects to be passed to other functions of the `gips`
package.

## Usage

``` r
gips_perm(x, size)

new_gips_perm(rearranged_cycles, size)

validate_gips_perm(g)
```

## Arguments

- x:

  An object created with a `permutations` package or any object that can
  be processed with the
  [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html)
  function.

- size:

  An integer. Size of a permutation (AKA cardinality of a set, on which
  permutation is defined; see examples).

- rearranged_cycles:

  A list of rearranged integer vectors. Each vector corresponds to a
  single cycle of a permutation.

- g:

  Object to be checked whether it is a proper object of a `gips_perm`
  class.

## Value

`gips_perm()` returns an object of a `gips_perm` class after the safety
checks.

`new_gips_perm()` returns an object of a `gips_perm` class without the
safety checks.

`validate_gips_perm()` returns its argument unchanged. If the argument
is not a proper element of a `gips_perm` class, it produces an error.

## Functions

- `new_gips_perm()`: Constructor. Only intended for low-level use.

- `validate_gips_perm()`: Validator. Only intended for low-level use.

## Methods for a `gips` class

- [`as.character.gips_perm()`](https://przechoj.github.io/gips/reference/as.character.md)

- [`print.gips_perm()`](https://przechoj.github.io/gips/reference/print.gips_perm.md)

## See also

- [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html) -
  The constructor for the `x` parameter.

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md) - The
  constructor for the `gips` class uses the `gips_perm` object as the
  base object.

## Examples

``` r
gperm <- gips_perm(permutations::as.word(c(1, 2, 3, 5, 4)), 5)
gperm <- gips_perm(permutations::as.cycle("(5,4)"), 5)
# note the necessity of `size` parameter
gperm <- gips_perm(permutations::as.cycle("(5,4)"), 7)
gperm <- gips_perm("(1,2)(5,4)", 7)
gperm
#> [1] (12)(45)

# \donttest{
try(gperm <- gips_perm(permutations::as.cycle("(5,4)"), 3))
#> Error in wrong_argument_abort(i = "`size` attribute must be greater or equal to largest integer in elements of `x`.",  : 
#>   There was a problem identified with provided argument
#> ℹ `size` attribute must be greater or equal to largest integer in elements of `x`.
#> ✖ `size` equals 3 while the maximum element is 5
# Error, `size` equals 3 while the maximum element is 5.
# }
```
