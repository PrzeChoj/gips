# Plot optimized matrix or optimization `gips` object

Plot the heatmap of the MAP covariance matrix estimator or the
convergence of the optimization method. The plot depends on the `type`
argument.

## Usage

``` r
# S3 method for class 'gips'
plot(
  x,
  type = NA,
  logarithmic_y = TRUE,
  logarithmic_x = FALSE,
  color = NULL,
  title_text = "Convergence plot",
  xlabel = NULL,
  ylabel = NULL,
  show_legend = TRUE,
  ylim = NULL,
  xlim = NULL,
  ...
)
```

## Arguments

- x:

  Object of a `gips` class.

- type:

  A single character. One of
  `c("heatmap", "block_heatmap", "all", "best", "both")`.

  - "heatmap" - Plots a heatmap of the `S` matrix inside the `gips`
    object projected on the permutation in the `gips` object.

  - "block_heatmap" - Plots a heatmap of diagonally block representation
    of `S`. Non-block entries (equal to 0) are white for better clarity.
    For more information see **Block Decomposition - \[1\], Theorem 1**
    section in
    [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.md)
    or in its [pkgdown
    page](https://przechoj.github.io/gips/articles/Theory.html)).

  - "all" - Plots the line of a posteriori for all visited states.

  - "best" - Plots the line of the biggest a posteriori found over time.

  - "both" - Plots both lines from "all" and "best".

  The default value is `NA`, which will be changed to "heatmap" for
  non-optimized `gips` objects and to "both" for optimized ones. Using
  the default produces a warning. All other arguments are ignored for
  the `type = "heatmap"`.

- logarithmic_y, logarithmic_x:

  A boolean. Sets the axis of the plot in logarithmic scale.

- color:

  Vector of colors to be used to plot lines.

- title_text:

  Text to be in the title of the plot.

- xlabel:

  Text to be on the bottom of the plot.

- ylabel:

  Text to be on the left of the plot.

- show_legend:

  A boolean. Whether or not to show a legend.

- ylim:

  Limits of the y axis. When `NULL`, the minimum and maximum of the
  [`log_posteriori_of_gips()`](https://przechoj.github.io/gips/reference/log_posteriori_of_gips.md)
  are taken.

- xlim:

  Limits of the x axis. When `NULL`, the whole optimization process is
  shown.

- ...:

  Additional arguments passed to
  [`stats::heatmap()`](https://rdrr.io/r/stats/heatmap.html) or other
  various elements of the plot.

## Value

Returns an invisible `NULL`.

## See also

- [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md) -
  Usually, the `plot.gips()` is called on the output of
  [`find_MAP()`](https://przechoj.github.io/gips/reference/find_MAP.md).

- [`project_matrix()`](https://przechoj.github.io/gips/reference/project_matrix.md) -
  The function used with `type = "heatmap"`.

- [`gips()`](https://przechoj.github.io/gips/reference/gips.md) - The
  constructor of a `gips` class. The `gips` object is used as the `x`
  parameter.

## Examples

``` r
require("MASS") # for mvrnorm()

perm_size <- 6
mu <- runif(6, -10, 10) # Assume we don't know the mean
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
number_of_observations <- 13
Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
S <- cov(Z) # Assume we have to estimate the mean

g <- gips(S, number_of_observations)
if (require("graphics")) {
  plot(g, type = "heatmap")
}


g_map <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "hill_climbing")
if (require("graphics")) {
  plot(g_map, type = "both", logarithmic_x = TRUE)
}


if (require("graphics")) {
  plot(g_map, type = "heatmap")
}

# Now, the output is (most likely) different because the permutation
  # `g_map[[1]]` is (most likely) not an identity permutation.
```
