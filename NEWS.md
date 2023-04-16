# gips 1.0.1.9052

### Bugfixes:

* Sometimes `post_probabilities` underflows to 0. This is appropriately validated now.
* `NaN`s should not occur in `find_MAP()` for `D_matrix <- diag(ncol(S)) * D_coef` when `1000 < D_coef < 1e300`.
* When `NaN`s do occur in `find_MAP()`, they will throw an error (used to show a warning).
* `Inf` better handled in `print.gips()`
* `print.*()` will print `\n` in the end
* Small Vignettes improvements
* Proper testing of examples


# gips 1.0.0

* This is the first release of gips.
