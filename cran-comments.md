## Submission

This is a maintenance update from gips 1.3.0 to 1.3.1. Changes are described
in NEWS.md.

## Test environment

* macOS 27.2 (ARM64), R 4.3.1.

## R CMD check results

Checked gips 1.3.1 with
`devtools::check(remote = TRUE, manual = TRUE)`
on macOS 27.2 (ARM64), R 4.3.1.

0 errors | 0 warnings | 0 notes

## GitHub Actions checks

The R CMD check workflow for gips 1.3.1 completed
successfully in all five environments:

* macOS-latest, R release.
* windows-latest, R release.
* ubuntu-latest, R release, R devel and R oldrel-1.

All five environments reported:
0 errors | 0 warnings | 0 notes

Results: <https://github.com/PrzeChoj/gips/actions/runs/42517049558>

## Minimum R version

The minimum R version declared in DESCRIPTION has been corrected to 4.0.0
to match existing code requirements, including class-specific warning
suppression.

gips 1.3.1 passed full R CMD build and R CMD check
runs on R 4.0.0 with C++11 and R 4.3.1 with C++17 on macOS ARM64,
using compatible dependencies, including archived versions where
needed for R 4.0.0. Both checks reported 0 errors, 0 warnings and
0 notes and included tests, examples, vignette rebuilding and the
PDF manual. These checks were run without --as-cran.

## Reverse dependencies

We checked gipsDA 1.0.0 on macOS 27.2 (ARM64), R 4.3.1,
comparing the CRAN version of gips (1.3.0) with the release candidate
(1.3.1). Both checks of gipsDA completed with Status: OK.

* Checked 1 reverse dependency.
* Found 0 new problems.
* Failed to check 0 packages.
