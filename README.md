
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gips`

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/gips)](https://CRAN.R-project.org/package=gips)
[![R-CMD-check](https://github.com/PrzeChoj/gips/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PrzeChoj/gips/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://codecov.io/gh/PrzeChoj/gips/branch/main/graph/badge.svg)](https://codecov.io/gh/PrzeChoj/gips?branch=main)
<!-- badges: end -->

gips - Gaussian model Invariant by Permutation Symmetry

`gips` is an R package that performs the Bayesian model selection
procedure within Gaussian models invariant by permutation symmetry
described by a cyclic subgroup. Permutation invariance of the Gaussian
distribution results in certain symmetries of its covariance matrix.
Such symmetry conditions reduce the number of parameters to estimate.
This is especially useful when parsimony is needed, i.e. when the number
of variables is substantially larger than the number of observations.
Given the Gaussian multivariate sample and two hyperparameters, through
a brute-force search or the Metropolis-Hasting algorithm, `gips` will
try to find a cyclic subgroup that best fits the data and will return
the estimated posterior probabilities for all cyclic subgroups.

## `gips` will help you with two things:

1.  Finding hidden symmetries between the variables. `gips` can be used
    as an exploratory tool searching the space of permutation symmetries
    of the Gaussian vector.
2.  Covariance estimation. The MLE for covariance matrix is known to
    exist if and only if the number of variables is less or equal to the
    number of observations. Additional knowledge of symmetries allows to
    significantly weaken this requirement. Moreover, the reduction of
    model dimension brings the advantage in terms of precision of
    covariance estimation.

## Installation

You can install the development version of gips from
[GitHub](https://github.com/PrzeChoj/gips) with:

``` r
# install.packages("devtools")
devtools::install_github("PrzeChoj/gips")
```

## Examples

### Example 1 - EDA

Assume we have the data, and we want to understand its structure:

``` r
library(gips)

Z <- as.matrix(mtcars)

# Assume the data is normal.
  # Looking at this (`hist(Z[,2])`) distribution,
  # it is not a remarkably sensible assumption,
  # but let's do it for the example.

S <- cor(Z)
g <- gips(S, nrow(Z), was_mean_estimated = TRUE)
plot(g, type = 'heatmap')
```

<img src="man/figures/README-example_mean_unknown-1.png" width="100%" />

``` r
# We can see some strong relationships between columns in this matrix.
  # For example, 9 and 10 have very similar correlations to other variables.

# Let's see if the find_MAP will find this relationship:
g_MAP <- find_MAP(g, max_iter = 10, optimizer = "MH")
#> ========================================================================
plot(g_MAP, type = 'heatmap')
```

<img src="man/figures/README-example_mean_unknown-2.png" width="100%" />

``` r
# Even after a short time (only 10 iterations),
  # find_MAP found some relationship.

# Let's see what it will find with a slightly bigger budget:
g_MAP <- find_MAP(g_MAP, max_iter = 100, optimizer = "continue")
#> ===============================================================================
plot(g_MAP, type = 'heatmap')
```

<img src="man/figures/README-example_mean_unknown-3.png" width="100%" />

``` r
# find_MAP found the (9,10) relationship and even something more.
```

### Example 2 - modeling

Assume we know the mean is 0, and we want to estimate the covariance
matrix, but we don’t have enough data:

``` r
# Prepare model, multivariate normal distribution
perm_size <- 6
mu <- numeric(perm_size)  
sigma_matrix <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)

# generate example data from a model:
Z <- MASS::mvrnorm(4, mu = mu, Sigma = sigma_matrix)
# End of prepare model



library(gips)
(number_of_observations <- nrow(Z)) # 4 < 6, so n < p
#> [1] 4

# calculate the covariance matrix from the data:
S <- (t(Z) %*% Z) / number_of_observations

# Make the gips object out of data:
g <- gips(S, number_of_observations, was_mean_estimated = FALSE)

# Find the Maximum A Posteriori Estimator for the permutation.
  # Space is small (6! = 720), so it is reasonable to
  # browse the whole of it:
g_map <- find_MAP(g, show_progress_bar = TRUE, optimizer = "full")
#> ================================================================================
g_map
#> The permutation (2,5,3,6) was found after 720 log_posteriori calculations, is 123.845008471236 times more likely than the starting, () permutation.

summary(g_map)$n0
#> [1] 3
summary(g_map)$n0 <= 4
#> [1] TRUE

# We see the number of observations (4) is bigger or equal to n0,
  # so we can estimate the covariance matrix
  # with the Maximum Likelihood estimator:
project_matrix(S, g_map[[1]])
#>              [,1]       [,2]       [,3]         [,4]       [,5]       [,6]
#> [1,]  0.098633812 0.04607666 0.04607666 -0.006480499 0.04607666 0.04607666
#> [2,]  0.046076656 1.16676492 1.08584461  1.352890482 0.27266235 0.27266235
#> [3,]  0.046076656 1.08584461 1.16676492  1.352890482 0.27266235 0.27266235
#> [4,] -0.006480499 1.35289048 1.35289048  2.712261501 1.35289048 1.35289048
#> [5,]  0.046076656 0.27266235 0.27266235  1.352890482 1.16676492 1.08584461
#> [6,]  0.046076656 0.27266235 0.27266235  1.352890482 1.08584461 1.16676492

# Plot the found matrix:
plot(g_map, type = "heatmap")
```

<img src="man/figures/README-example_mean_known-1.png" width="100%" />

# Credits

`gips` was developed in 2022 by Przemysław Chojecki and Paweł Morgen
under the leadership of Ph.D. Bartosz Kołodziejek within the
“CyberiADa-3” (2021) grant from the Warsaw University of Technology.
