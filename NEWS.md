# gips 1.0.1.9000

### Bugfixes:

* Sometimes `post_probabilities` underflows to 0. This is appropriately validated now.
* Typos in Vignettes
* Proper testing of examples
* When `NaN`s occur in `find_MAP`, they will throw an error (used to show a warning).
* `NaN`s should not occur for `D_matrix <- diag(ncol(S)) * D_coef` when `1000 < D_coef < 1e300`.


# gips 1.0.0

* This is the first release of gips.
