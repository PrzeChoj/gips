# gips 1.1.0.9103

### New functions

-   `BIC.gips()`
-   `AIC.gips()`
-   `logLik.gips()`
-   `as.character.gips()`

### Update to functions

-   `gips()` has a new default `D_matrix = mean(diag(S)) * diag(p)`;
-   `summary.gips()` calculates `AIC`, `BIC`, and `n_parameters` (number of free parameters in the covariance matrix);
-   `get_probabilities_from_gips()` will return a sorted vector;
-   `compare_posteriories_of_perms()` and `compare_log_posteriories_of_perms()` have a new parameter `digits`;
-   Everywhere a permutation was expected, the `gips` object can now be passed and interpreted as a permutation. Those are:
    -   `perm` in `gips()`, `project_matrix()`, `prepare_orthogonal_matrix()`, `get_structure_constants()`, `calculate_gamma_function()`;
    -   `perm1` and `perm2` in `compare_posteriories_of_perms()`, `compare_log_posteriories_of_perms()`;
    -   `x` in `gips_perm()`;
-   `plot.gips()` can get `type = "MLE"`, which is an alias for `type = "heatmap"`;
-   `find_MAP(optimizer = "BF")` is 3 times faster;
-   `find_MAP(optimizer = "BF")` is default for `perm_size <= 9`.

### Bugfixes:

-   Sometimes, `post_probabilities` underflows to 0. This is appropriately validated now;
-   `NaN`s should not occur in `find_MAP()` for `D_matrix <- diag(ncol(S)) * d` when `1000 < d < 1e300`;
-   When `NaN`s do occur in `find_MAP()`, they will throw an error (used to show a warning);
-   `Inf` better handled in `print.gips()`;
-   `print.*()` functions will print `\n` in the end;
-   Slightly different punctuation in `print.gips()`;
-   Tremendous vignettes and documentation improvements;
-   Proper testing of examples;
-   `delta` parameter of `gips()` has to be bigger than `1`. We used to restrict it to bigger than `2`;
-   `project_matrix()` shows a warning when the non-positive-semi-definite matrix is passed as an `S` argument;
-   `project_matrix()` preserves `colnames()` and `rownames()` of a matrix;
-   `D_matrix` is checked for containing any `NaN` or `Inf` values;
-   Absurdly long structure constants vectors may overflow an `integer`. Now we use `double`;
-   `compare_log_posteriories_of_perms()` and `compare_posteriories_of_perms()` show an error when given two incomparable `gips` objects (with different parameters).

# gips 1.0.0

-   This is the first release of gips.
