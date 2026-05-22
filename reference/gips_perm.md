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

  A single object that can be interpreted by the
  [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html)
  function. For example, the character of a form `"(1,2)(4,5)"`. See
  examples. It can also be of a `gips` class but it will be interpreted
  as the underlying `gips_perm`.

- size:

  An integer. Size of a permutation (AKA cardinality of a set, on which
  permutation is defined. See examples).

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

- [`as.character.gips_perm()`](https://przechoj.github.io/gips/reference/as.character.gips_perm.md)

- [`print.gips_perm()`](https://przechoj.github.io/gips/reference/print.gips_perm.md)

## See also

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  `gips_perm` is the `perm` parameter of
  [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md).

- [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html) -
  The constructor for the `x` parameter.

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md) - The
  constructor for the `gips` class uses the `gips_perm` object as the
  base object.

## Examples

``` r
# All 7 following lines give the same output:
gperm <- gips_perm("(12)(45)", 5)
gperm <- gips_perm("(1,2)(4,5)", 5)
gperm <- gips_perm(as.matrix(c(2, 1, 3, 5, 4)), 5)
gperm <- gips_perm(t(as.matrix(c(2, 1, 3, 5, 4))), 5) # both way for a matrix works
gperm <- gips_perm(list(list(c(2, 1), c(4, 5))), 5)
gperm <- gips_perm(permutations::as.word(c(2, 1, 3, 5, 4)), 5)
gperm <- gips_perm(permutations::as.cycle("(1,2)(4,5)"), 5)
gperm
#> [1] (12)(45)

# note the necessity of the `size` parameter:
gperm <- gips_perm("(12)(45)", 5)
gperm <- gips_perm("(12)(45)", 7) # this one is a different permutation

try(gperm <- gips_perm("(12)(45)", 4))
#> Error in wrong_argument_abort(i = "`size` attribute must be greater or equal to the largest integer in elements of `x`.",  : 
#>   There was a problem identified with provided argument
#> ℹ `size` attribute must be greater or equal to the largest integer in elements of `x`.
#> ✖ `size` equals 4 while the maximum element is 5
# Error, `size` was set to 4, while the permutation has the element 5.
```
