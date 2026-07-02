# Find the Maximum A Posteriori Estimation

Use one of the optimization algorithms to find the permutation that
maximizes A Posteriori probability based on observed data. Not all
optimization algorithms will always find the MAP, but they try to find a
significant value. More information can be found in the "**Possible
algorithms to use as optimizers**" section below.

## Usage

``` r
find_MAP(
  g,
  max_iter = NA,
  optimizer = NA,
  show_progress_bar = TRUE,
  save_all_perms = FALSE,
  return_probabilities = FALSE
)
```

## Arguments

- g:

  Object of a `gips` class.

- max_iter:

  The number of iterations for an algorithm to perform. At least 2. For
  `optimizer = "BF"`, it is not used; for `optimizer = "MH"`, it has to
  be finite; for `optimizer = "HC"`, it can be infinite.

- optimizer:

  The optimizer for the search of the maximum posteriori:

  - `"BF"` (the default for unoptimized `g` with `perm size <= 9`) -
    Brute Force;

  - `"MH"` (the default for unoptimized `g` with `perm size > 10`) -
    Metropolis-Hastings;

  - `"HC"` - Hill Climbing;

  - `"continue"` (the default for optimized `g`) - The same as the `g`
    was optimized by (see Examples).

  See the **Possible algorithms to use as optimizers** section below for
  more details.

- show_progress_bar:

  A boolean. Indicate whether or not to show the progress bar:

  - When `max_iter` is infinite, `show_progress_bar` has to be `FALSE`;

  - When `return_probabilities = TRUE`, then shows an additional
    progress bar for the time when the probabilities are calculated.

- save_all_perms:

  A boolean. `TRUE` indicates saving a list of all permutations visited
  during optimization. This can be useful sometimes but needs a lot more
  RAM.

- return_probabilities:

  A boolean. `TRUE` can only be provided when `save_all_perms = TRUE`.
  For:

  - `optimizer = "MH"` - use Metropolis-Hastings results to estimate
    posterior probabilities;

  - `optimizer = "BF"` - use brute force results to calculate exact
    posterior probabilities.

  These additional calculations are costly, so a second and third
  progress bars are shown (when `show_progress_bar = TRUE`).

  To examine probabilities after optimization, call
  [`get_probabilities_from_gips()`](https://przechoj.github.io/gips/reference/get_probabilities_from_gips.md).

## Value

An optimized object of a `gips` class.

## Details

`find_MAP()` can produce a warning when:

- the optimizer "hill_climbing" gets to the end of its `max_iter`
  without converging.

- the optimizer will find the permutation with smaller `n0` than
  `number_of_observations` (for more information on what it means, see
  **\\C\_\sigma\\ and `n0`** section in the
  [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
  or in its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html)).

## Possible algorithms to use as optimizers

For an in-depth explanation, see in the
[`vignette("Optimizers", package = "gips")`](https://przechoj.github.io/gips/articles/Optimizers.md)
or in its [pkgdown
page](https://przechoj.github.io/gips/articles/Optimizers.html).

For every algorithm, there are some aliases available.

- `"brute_force"`, `"BF"`, `"full"` - use the **Brute Force** algorithm
  that checks the whole permutation space of a given size. This
  algorithm will find the actual Maximum A Posteriori Estimation, but it
  is very computationally expensive for bigger spaces. We recommend
  Brute Force only for `p <= 9`. For the time the Brute Force takes on
  our machines, see in the
  [`vignette("Optimizers", package = "gips")`](https://przechoj.github.io/gips/articles/Optimizers.md)
  or in its [pkgdown
  page](https://przechoj.github.io/gips/articles/Optimizers.html).

- `"Metropolis_Hastings"`, `"MH"` - use the **Metropolis-Hastings**
  algorithm; [see
  Wikipedia](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
  The algorithm will draw a random transposition in every iteration and
  consider changing the current state (permutation). When the `max_iter`
  is reached, the algorithm will return the best permutation calculated
  as the MAP Estimator. This implements the [*Second approach* from
  references, section 4.1.2](https://arxiv.org/abs/2004.03503). This
  algorithm used in this context is a special case of the **Simulated
  Annealing** the user may be more familiar with; [see
  Wikipedia](https://en.wikipedia.org/wiki/Simulated_annealing).

- `"hill_climbing"`, `"HC"` - use the **hill climbing** algorithm; [see
  Wikipedia](https://en.wikipedia.org/wiki/Hill_climbing). The algorithm
  will check all transpositions in every iteration and go to the one
  with the biggest A Posteriori value. The optimization ends when all
  *neighbors* will have a smaller A Posteriori value. If the `max_iter`
  is reached before the end, then the warning is shown, and it is
  recommended to continue the optimization on the output of the
  `find_MAP()` with `optimizer = "continue"`; see examples. Remember
  that `p*(p-1)/2` transpositions will be checked in every iteration.
  For bigger `p`, this may be costly.

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md) - The
  constructor of a `gips` class. The `gips` object is used as the `g`
  parameter of `find_MAP()`.

- [`plot.gips()`](https://przechoj.github.io/gips/reference/plot.gips.md) -
  Practical plotting function for visualizing the optimization process.

- [`summary.gips()`](https://przechoj.github.io/gips/reference/summary.gips.md) -
  Summarize the output of optimization.

- [`AIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md),
  [`BIC.gips()`](https://przechoj.github.io/gips/reference/AIC.gips.md) -
  Get the Information Criterion of the found model.

- [`get_probabilities_from_gips()`](https://przechoj.github.io/gips/reference/get_probabilities_from_gips.md) -
  When `find_MAP(return_probabilities = TRUE)` was called, probabilities
  can be extracted with this function.

- [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md) -
  The function that the optimizers of `find_MAP()` tries to find the
  argmax of.

- [`forget_perms()`](https://przechoj.github.io/gips/reference/forget_perms.md) -
  When the `gips` object was optimized with
  `find_MAP(save_all_perms = TRUE)`, it will be of considerable size in
  RAM.
  [`forget_perms()`](https://przechoj.github.io/gips/reference/forget_perms.md)
  can make such an object lighter in memory by forgetting the
  permutations it visited.

- [`vignette("Optimizers", package = "gips")`](https://przechoj.github.io/gips/articles/Optimizers.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Optimizers.html) - A
  place to learn more about the available optimizers.

- [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html) - A place
  to learn more about the math behind the `gips` package.

## Examples

``` r
require("MASS") # for mvrnorm()

perm_size <- 5
mu <- runif(perm_size, -10, 10) # Assume we don't know the mean
sigma_matrix <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5)
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)

g_map <- find_MAP(g, max_iter = 5, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
g_map
#> The permutation (1,2):
#>  - was found after 5 posteriori calculations;
#>  - is 2.337 times more likely than the () permutation.

g_map2 <- find_MAP(g_map, max_iter = 5, show_progress_bar = FALSE, optimizer = "continue")
plot(g_map2, type = "both", logarithmic_x = TRUE)


g_map_BF <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
summary(g_map_BF)
#> The optimized `gips` object.
#> 
#> Permutation:
#>  (1,2,3,4,5)
#> 
#> Log_posteriori:
#>  -15.1521
#> 
#> Times more likely than starting permutation:
#>  414.386
#> 
#> The p-value of Likelihood-Ratio test:
#>  0.4883
#> 
#> The number of observations:
#>  13
#> 
#> The mean in the `S` matrix was estimated.
#> Therefore, one degree of freedom was lost.
#> There are 12 degrees of freedom left.
#> 
#> n0:
#>  2
#> 
#> The number of observations is bigger than n0 for this permutation,
#> so the gips model based on the found permutation does exist.
#> 
#> The number of free parameters in the covariance matrix:
#>  3
#> 
#> BIC:
#>  129.0236
#> 
#> AIC:
#>  127.3288
#> 
#> --------------------------------------------------------------------------------
#> Optimization algorithm:
#>  brute_force
#> 
#> The number of log_posteriori calls:
#>  67
#> 
#> Optimization time:
#>  0.07670975 secs
```
