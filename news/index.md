# Changelog

## gips 1.2.3.9000

#### Update to functions

- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md)
  now uses `ggplot2` for all plot types.

#### Bugfix:

- Documentation improvements: grammar and style corrections in roxygen
  comments, vignettes, and error messages.

## gips 1.2.3

CRAN release: 2025-03-18

#### Performance gain

There was a significant speed improvement in
[`get_structure_constants()`](https://przechoj.github.io/gips/reference/get_structure_constants.md),
which is used internally by posterior calculations. For 1000 random
permutations of a given size:

| permutation size |     30 |     50 |    100 |    200 |     300 |
|------------------|-------:|-------:|-------:|-------:|--------:|
| v1.2.2           | 0.09 s | 0.10 s | 0.25 s |  ~10 s |   ~25 s |
| v1.2.3           | 0.07 s | 0.08 s | 0.10 s | 0.17 s | ~0.20 s |

This change consequently improved
[`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md).
For 1000 random permutations of a given size:

| permutation size |     30 |     50 |    100 |     200 |      300 |
|------------------|-------:|-------:|-------:|--------:|---------:|
| v1.2.2           | 1.50 s | 2.49 s | 7.08 s | 55.25 s | 169.77 s |
| v1.2.3           | 1.29 s | 2.24 s | 6.84 s | 31.73 s |  85.10 s |

#### Update to functions

- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md)
  can get `type = "n0"`, which will plot the change of `n0` along the
  “MH” optimization. Handy for deciding of burn-in time;
- `find_MAP(optimizer = "MH")` tracks the `n0` along the optimization;
- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md)
  calculates Likelihood-Ratio test.

## gips 1.2.2

#### Bugfix:

- [`logLik.gips()`](https://przechoj.github.io/gips/reference/logLik.gips.md)
  will return an object of class `logLik`;
- Better Vignettes titles.

## gips 1.2.1

CRAN release: 2023-08-12

#### Bugfix:

- One test on alternative BLAS/LAPACK has an alternative check for a
  singular matrix. Found by CRAN’s ATLAS test.

## gips 1.2.0

CRAN release: 2023-08-07

#### New functions

- [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
- [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md)
- [`logLik.gips()`](https://przechoj.github.io/gips/reference/logLik.gips.md)
- [`as.character.gips()`](https://przechoj.github.io/gips/reference/as.character.gips.md)

#### Update to functions

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md) has a
  new default `D_matrix = mean(diag(S)) * diag(p)`;
- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md)
  calculates `AIC`, `BIC`, and `n_parameters` (number of free parameters
  in the covariance matrix);
- [`get_probabilities_from_gips()`](https://przechoj.github.io/gips/reference/get_probabilities_from_gips.md)
  will return a sorted vector;
- [`compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  and
  [`compare_log_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  have a new parameter `digits`;
- Everywhere a permutation was expected, the `gips` object can now be
  passed and interpreted as a permutation. Those are:
  - `perm` in
    [`gips()`](https://przechoj.github.io/gips/reference/gips.md),
    [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md),
    [`prepare_orthogonal_matrix()`](https://przechoj.github.io/gips/reference/prepare_orthogonal_matrix.md),
    [`get_structure_constants()`](https://przechoj.github.io/gips/reference/get_structure_constants.md),
    [`calculate_gamma_function()`](https://przechoj.github.io/gips/reference/calculate_gamma_function.md);
  - `perm1` and `perm2` in
    [`compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md),
    [`compare_log_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md);
  - `x` in
    [`gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.md);
- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md)
  can get `type = "MLE"`, which is an alias for `type = "heatmap"`;
- `find_MAP(optimizer = "BF")` is 3 times faster;
- `find_MAP(optimizer = "BF")` is default for `perm_size <= 9`.

#### Bugfixes:

- Sometimes, `post_probabilities` underflows to 0. This is appropriately
  validated now;
- `NaN`s should not occur in
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md)
  for `D_matrix <- diag(ncol(S)) * d` when `1000 < d < 1e300`;
- When `NaN`s do occur in
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md),
  they will throw an error (used to show a warning);
- `Inf` better handled in
  [`print.gips()`](https://przechoj.github.io/gips/reference/print.gips.md);
- `print.*()` functions will print `\n` in the end;
- Slightly different punctuation in
  [`print.gips()`](https://przechoj.github.io/gips/reference/print.gips.md);
- Tremendous vignettes and documentation improvements;
- Proper testing of examples;
- `delta` parameter of
  [`gips()`](https://przechoj.github.io/gips/reference/gips.md) has to
  be bigger than `1`. We used to restrict it to bigger than `2`;
- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
  shows a warning when the non-positive-semi-definite matrix is passed
  as an `S` argument;
- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md)
  preserves [`colnames()`](https://rdrr.io/r/base/colnames.html) and
  [`rownames()`](https://rdrr.io/r/base/colnames.html) of a matrix;
- `D_matrix` is checked for containing any `NaN` or `Inf` values;
- Absurdly long structure constants vectors may overflow an `integer`.
  Now we use `double`;
- [`compare_log_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  and
  [`compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.md)
  show an error when given two incomparable `gips` objects (with
  different parameters).

## gips 1.0.0

CRAN release: 2022-10-13

- This is the first release of gips.
