# Package index

## Main functions

Functions for the main `gips` purpose: Finding the permutation under
which a matrix is (approximately) invariant and estimating the
covariance matrix.

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md)
  [`new_gips()`](https://przechoj.github.io/gips/reference/gips.md)
  [`validate_gips()`](https://przechoj.github.io/gips/reference/gips.md)
  :

  The constructor of a `gips` class.

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
  : Find the Maximum A Posteriori Estimation

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
  : Project matrix after optimization

## Diagnostics

Functions for diagnostics of the `gips` optimization process or its
results.

- [`plot(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/plot.gips.md)
  :

  Plot optimized matrix or optimization `gips` object

- [`print(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/print.gips.md)
  :

  Printing `gips` object

- [`summary(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/summary.gips.md)
  [`print(`*`<summary.gips>`*`)`](https://przechoj.github.io/gips/reference/summary.gips.md)
  : Summarizing the gips object

- [`compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  [`compare_log_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  : Compare the posteriori probabilities of 2 permutations

- [`AIC(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/AIC.gips.md)
  [`BIC(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/AIC.gips.md)
  :

  Akaike's An Information Criterion for `gips` class

- [`logLik(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/logLik.gips.md)
  :

  Extract the Log-Likelihood for `gips` class

## Helper functions

Functions for pleasant work in `gips`.

- [`as.character(`*`<gips>`*`)`](https://przechoj.github.io/gips/reference/as.character.gips.md)
  :

  Transform the `gips` object to a character vector

- [`gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md)
  [`new_gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md)
  [`validate_gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md)
  : Permutation object

- [`print(`*`<gips_perm>`*`)`](https://przechoj.github.io/gips/reference/print.gips_perm.md)
  :

  Printing `gips_perm` object

- [`as.character(`*`<gips_perm>`*`)`](https://przechoj.github.io/gips/reference/as.character.gips_perm.md)
  :

  Transform the `gips_perm` object to a character vector

- [`forget_perms()`](https://przechoj.github.io/gips/reference/forget_perms.md)
  :

  Forget the permutations for `gips` object optimized with
  `save_all_perms = TRUE`

- [`get_probabilities_from_gips()`](https://przechoj.github.io/gips/reference/get_probabilities_from_gips.md)
  :

  Extract probabilities for `gips` object optimized with
  `return_probabilities = TRUE`

## Internal functions

Functions for internal calculations

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
  : Log of the A Posteriori of the permutation that the covariance
  matrix is invariant under
- [`prepare_orthogonal_matrix()`](https://przechoj.github.io/gips/reference/prepare_orthogonal_matrix.md)
  : Prepare orthogonal matrix
- [`get_structure_constants()`](https://przechoj.github.io/gips/reference/get_structure_constants.md)
  : Get Structure Constants
- [`calculate_gamma_function()`](https://przechoj.github.io/gips/reference/calculate_gamma_function.md)
  : Calculate Gamma function
