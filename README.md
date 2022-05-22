
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gips`

<!-- badges: start -->
<!-- badges: end -->

gips - Gaussian model Invariant by Permutation Symmetry

`gips` is an R package that finds the permutation symmetry group such
that the covariance matrix of the given data is invariant under it.
Knowledge of such a permutation can drastically decrease the number of
parameters needed to fit the model. That means that with `gips`, it is
possible to find the Gaussian model with more parameters than the number
of observations. Sometimes, even if the number of observations is bigger
than the number of parameters, the covariance matrix found with `gips`
better approximates the actual covariance behind the data.

## `gips` will help you with two things:

1.  Exploratory Data Analysis - with `gips`, you can find the
    permutation of features that does not change the covariance matrix.
2.  Modeling - with `gips`, you can accurately use the found permutation
    to fit the normal models like LDA or QDA.

## Installation

You can install the development version of gips from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("PrzeChoj/gips")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(gips)
## TODO(too little data to estimate covariance matrix -> `gips` will find the matrix)
```

# Credits

It is developed by Przemysław Chojecki and Paweł Morgen under the
leadership of Ph.D. Bartosz Kołodziejek.
