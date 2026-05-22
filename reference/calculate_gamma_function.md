# Calculate Gamma function

It calculates the value of the integral defined in [Definition 11 from
references](https://arxiv.org/abs/2004.03503). It implements [Theorem 8
from references](https://arxiv.org/abs/2004.03503) and uses the [formula
(19) from references](https://arxiv.org/abs/2004.03503).

## Usage

``` r
calculate_gamma_function(perm, lambda)
```

## Arguments

- perm:

  An object of a `gips_perm` class. It can also be of a `gips` class,
  but it will be interpreted as the underlying `gips_perm`.

- lambda:

  A positive real number.

## Value

Returns the value of the Gamma function of the colored cone (for the
definition of the colored cone, see the **Basic definitions** section in
[`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or in its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html)).

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- [`get_structure_constants()`](https://przechoj.github.io/gips/reference/get_structure_constants.md) -
  The function useful inside the `calculate_gamma_function()`.

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  The function that uses the values of the gamma function.

- [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html) - A place
  to learn more about the math behind the `gips` package.

## Examples

``` r
id_perm <- gips_perm("()", 2)
calculate_gamma_function(id_perm, 0.5001) # 10.7...
#> [1] 10.70139
calculate_gamma_function(id_perm, 0.50000001) # 19.9...
#> [1] 19.91198
calculate_gamma_function(id_perm, 0.500000000001) # 29.1...
#> [1] 29.12235

oldw <- getOption("warn")
options(warn = -1)
calculate_gamma_function(id_perm, 0.5) # Inf
#> [1] Inf
# Integral diverges; returns Inf and warning
options(warn = oldw)
```
