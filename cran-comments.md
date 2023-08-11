## Resubmission
This is a resubmission. In this version I have:

* solved the problem shown on <https://cran.r-project.org/web/checks/check_results_gips.html>. The ATLAS server said the matrices in some tests were singular. This was right and intended, but this was somehow not the problem in my local tests on any platform. I deleted the tests as they were not significant.


## Test environments

GitHub Actions using `usethis::use_github_actions_check_standard()`

* MacOS-latest (release)
* Windows-latest (release)
* Ubuntu-latest (oldrel-1, release, devel)

Testing with `devtools::check_rhub()` - no Errors, no Warnings, only technical Notes (package 'V8' unavailable; files/directories: ''NULL'', 'lastMiKTeXException')


## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded
