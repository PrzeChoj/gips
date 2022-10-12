## Resubmission
This is a resubmission. In this version I have changed the following:

* Used `\donttest{}` instead of `\dontrun{}` in man/gips_perm.Rd and man/calculate_gamma_function.Rd.
* Double spaces changed to single spaces in the description.


## Test environments

GitHub Actions using `usethis::use_github_actions_check_standard()`

* MacOS-latest (release)
* Windows-latest (release)
* Ubuntu-latest (oldrel-1, release, devel)

Testing with `devtools::check_rhub()` - no Errors, no Warnings, only technical Notes (Possibly misspelled words; URL doi; folder 'lastMiKTeXException')


## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

* This is the initial submission of this package to CRAN.
* This is my first submission of any package to CRAN.
* All exported functions document their return value.
* All exported functions documentations contain some runnable examples.
* We wrapped the examples of "print.*" functions in `\donttest{}`, so those will not be printed on the console while testing. We also wrapped the examples `man/gips_perm.Rd` and `man/calculate_gamma_function.Rd` in `\donttest{}`, because they purposefully shows errors after the improper use.
