## Resubmission
This is a resubmission. In this version I have:

* Change the URL in README.md from https://codecov.io/gh/PrzeChoj/gips?branch=main to https://app.codecov.io/gh/PrzeChoj/gips?branch=main

## Test environments

GitHub Actions using `usethis::use_github_actions_check_standard()`

* MacOS-latest (release)
* Windows-latest (release)
* Ubuntu-latest (oldrel-1, release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

* This is initial submission of this package to CRAN.
* This is my first submission of any package to CRAN.
* All exported functions document their return value.
* All exported functions documentations contain some runnable examples. We commented out examples of "print.*" functions so that those will not be printed on the console while testing.

## devtools::check_rhub()
Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1214/22-AOS2174
This is mysterious for me, because in the whole of documentation we used `\doi{10.1214/22-AOS2174}`.
