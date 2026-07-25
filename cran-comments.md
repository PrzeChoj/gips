## Test environments

GitHub Actions using the standard R CMD check workflow:

* macOS-latest (release)
* windows-latest (release)
* ubuntu-latest (oldrel-1)
* ubuntu-latest (release)
* ubuntu-latest (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resolved CRAN check issues

This release resolves both issues reported for the CRAN version of gips
(1.2.3):

* Replaced the deprecated `.Dim` argument to `structure()` with `dim`,
  resolving the R-devel NOTE.
* Made the exact-probability test independent of the relative ordering of tied
  probabilities, resolving the test failure under Intel MKL.

## Reverse dependencies

We checked the one reverse dependency, `gipsDA`, comparing R CMD check results
with the CRAN version of gips (1.2.3) and this version (1.3.0).

* We saw 0 new problems.
* We failed to check 0 packages.
