# gips 1.1.0.9000

### New functions

-   `BIC.gips()`
-   `AIC.gips()`
-   `logLik.gips()`

### Update to functions

-   `summary.gips()` calculates `AIC`, `BIC`, and `n_parameters` (number of free parameters in the covariance matrix).
-   `get_probabilities_from_gips()` has a new parameter `sorted`.
-   `compare_posteriories_of_perms()` and `compare_log_posteriories_of_perms()` have a new parameter `digits`.
-   Everywhere a permutation was expected, the `gips` object can now be passed and interpreted as a permutation. Those are:
    -   `perm` in `gips()`, `project_matrix()`, `prepare_orthogonal_matrix()`, `get_structure_constants()`, `calculate_gamma_function()`;
    -   `perm1` and `perm2` in `compare_posteriories_of_perms()`, `compare_log_posteriories_of_perms()`;
    -   `x` in `gips_perm()`.

### Bugfixes:

-   Sometimes, `post_probabilities` underflows to 0. This is appropriately validated now.
-   `NaN`s should not occur in `find_MAP()` for `D_matrix <- diag(ncol(S)) * d` when `1000 < d < 1e300`.
-   When `NaN`s do occur in `find_MAP()`, they will throw an error (used to show a warning).
-   `Inf` better handled in `print.gips()`
-   `print.*()` will print `\n` in the end
-   Small Vignettes and documentation improvements
-   Proper testing of examples
-   `delta` parameter of `gips()` has to be bigger than `1`. We used to restrict it as bigger than `2`.
-   `project_matrix()` will show a warning when the non-positive-semi-definite matrix is passed as an `S` argument.
-   `project_matrix()` will preserve colnames and rownames of a matrix.
-   `D_matrix` will be checked for containing any `NaN` or `Inf` values.
-   `plot.gips()` can get `type = "MLE"`, which is an alias for `type = "heatmap"`.


# gips 1.0.0

-   This is the first release of gips.
