# Get Structure Constants

Finds constants that are necessary for internal calculations of
integrals and eventually the posteriori probability in
[`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).

## Usage

``` r
get_structure_constants(perm)
```

## Arguments

- perm:

  An object of a `gips_perm` class.

## Value

Returns a list of 5 items: `r`, `d`, `k`, `L`, `dim_omega` - vectors of
constants from [Theorem 1 from
references](https://arxiv.org/abs/2004.03503) and the beginning of
[section 3.1. from references](https://arxiv.org/abs/2004.03503).

## Details

Uses the [Theorem 5 from references](https://arxiv.org/abs/2004.03503)
to calculate the constants.

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- [`calculate_gamma_function()`](https://przechoj.github.io/gips/reference/calculate_gamma_function.md),
  [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  The functions that rely heavily on `get_structure_constants()`.

## Examples

``` r
perm <- gips_perm(permutations::as.word(c(1, 2, 3, 5, 4)), 5)
get_structure_constants(perm)
#> $r
#> [1] 4 1
#> 
#> $d
#> [1] 1 1
#> 
#> $k
#> [1] 1 1
#> 
#> $L
#> [1] 2
#> 
#> $dim_omega
#> [1] 10  1
#> 
```
