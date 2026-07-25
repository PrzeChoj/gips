# Prepare orthogonal matrix

Calculate the orthogonal matrix `U_Gamma` for decomposition in [Theorem
1 from references](https://arxiv.org/abs/2004.03503).

## Usage

``` r
prepare_orthogonal_matrix(perm, perm_size = NULL, basis = NULL)
```

## Arguments

- perm:

  An object of a `gips_perm` or anything a
  [`gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md)
  can handle. It can also be of a `gips` class, but it will be
  interpreted as the underlying `gips_perm`.

- perm_size:

  Size of a permutation. Required if `perm` is neither `gips_perm` nor
  `gips`.

- basis:

  A matrix with basis vectors in COLUMNS. Identity by default.

## Value

A square matrix of size `perm_size` by `perm_size` with columns from
vector elements \\v_k^{(c)}\\ according to [Theorem 6 from
references](https://arxiv.org/abs/2004.03503).

## Details

Given X - a matrix invariant under the permutation `perm`. Call Gamma
the permutations cyclic group: \\\Gamma = \<perm\> = \\perm, perm^2,
...\\\\.

Then, \\U\_\Gamma\\ is such an orthogonal matrix, which
block-diagonalizes X.

To be more precise, the matrix `t(U_Gamma) %*% X %*% U_Gamma` has a
block-diagonal structure, which is ensured by [Theorem 1 from
references](https://arxiv.org/abs/2004.03503).

The formula for `U_Gamma` can be found in [Theorem 6 from
references](https://arxiv.org/abs/2004.03503).

A nice example is demonstrated in the **Block Decomposition - \[1\],
Theorem 1** section of
[`vignette("Theory", package="gips")`](https://przechoj.github.io/gips/articles/Theory.md)
or its [pkgdown
page](https://przechoj.github.io/gips/articles/Theory.html).

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  A function used in examples to show the properties of
  `prepare_orthogonal_matrix()`.

- **Block Decomposition - \[1\], Theorem 1** section of
  [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html) - A place
  to learn more about the math behind the `gips` package and see more
  examples of `prepare_orthogonal_matrix()`.

## Examples

``` r
gperm <- gips_perm("(1,2,3)(4,5)", 5)
U_Gamma <- prepare_orthogonal_matrix(gperm)

number_of_observations <- 10
X <- matrix(rnorm(5 * number_of_observations), number_of_observations, 5)
S <- cov(X)
X <- project_matrix(S, perm = gperm) # this matrix in invariant under gperm

block_decomposition <- t(U_Gamma) %*% X %*% U_Gamma
round(block_decomposition, 5) # the non-zeros only on diagonal and [1,2] and [2,1]
#>         [,1]    [,2]    [,3]    [,4]   [,5]
#> [1,] 0.41684 0.11768 0.00000 0.00000 0.0000
#> [2,] 0.11768 1.09027 0.00000 0.00000 0.0000
#> [3,] 0.00000 0.00000 0.92989 0.00000 0.0000
#> [4,] 0.00000 0.00000 0.00000 0.92989 0.0000
#> [5,] 0.00000 0.00000 0.00000 0.00000 0.9588
```
