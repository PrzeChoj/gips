#' Plot optimized matrix or optimization `gips` object
#'
#' Plot the heatmap of the MAP covariance matrix estimator
#' or the convergence of the optimization method.
#' The plot depends on the `type` argument.
#'
#' @param x Object of a `gips` class.
#' @param type A character vector of length 1. One of
#'     `c("heatmap", "MLE", "best", "all", "both", "n0", "block_heatmap")`:
#'   * `"heatmap"`, `"MLE"` - Plots a heatmap of the Maximum Likelihood
#'       Estimator of the covariance matrix given the permutation.
#'       That is, the `S` matrix inside the `gips` object
#'       projected on the permutation in the `gips` object.
#'   * `"best"` - Shows the maximum A Posteriori value found over time.
#'   * `"all"` - Shows the A Posteriori values for all visited states.
#'   * `"both"` - Shows both trajectories from "all" and "best".
#'   * `"n0"` - Plots the `n0` values observed during optimization
#'       (only for "MH" optimization).
#'   * `"block_heatmap"` - Plots a heatmap of diagonally block representation of `S`.
#'       Non-block entries (equal to 0) are white for better clarity.
#'       For more information, see **Block Decomposition - \[1\], Theorem 1**
#'       section in `vignette("Theory", package = "gips")` or in its
#'       [pkgdown page](https://przechoj.github.io/gips/articles/Theory.html).
#'
#' The default value is `NA`, which will be changed to "heatmap" for
#'     non-optimized `gips` objects and to "both" for optimized ones.
#'     Using the default produces a warning.
#'     All other arguments are ignored for
#'     the `type = "heatmap"`, `type = "MLE"`, or `type = "block_heatmap"`.
#' @param logarithmic_y,logarithmic_x A boolean.
#'     Sets the axis of the plot in logarithmic scale.
#' @param color Vector of colors to be used to plot lines.
#' @param title_text Text to be in the title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend A boolean. Whether or not to show a legend.
#' @param ylim Limits of the y axis. When `NULL`,
#'     the minimum, and maximum of the [log_posteriori_of_gips()] are taken.
#' @param xlim Limits of the x axis. When `NULL`,
#'     the whole optimization process is shown.
#' @param ... Additional arguments passed to
#'     other various elements of the plot.
#'
#' @returns When `type` is one of `"best"`, `"all"`, `"both"` or `"n0"`,
#'     returns an invisible `NULL`.
#'     When `type` is one of `"heatmap"`, `"MLE"` or `"block_heatmap"`,
#'     returns an object of class `ggplot`.
#'
#' @seealso
#' * [find_MAP()] - Usually, the `plot.gips()`
#'     is called on the output of `find_MAP()`.
#' * [project_matrix()] - The function used with `type = "MLE"`.
#' * [gips()] - The constructor of a `gips` class.
#'     The `gips` object is used as the `x` parameter.
#'
#' @export
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu <- runif(6, -10, 10) # Assume we don't know the mean
#' sigma_matrix <- matrix(
#'   data = c(
#'     1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'     0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'     0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'     0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'     0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'     0.8, 0.6, 0.4, 0.6, 0.8, 1.0
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' ) # sigma_matrix is a matrix invariant under permutation (1,2,3,4,5,6)
#' number_of_observations <- 13
#' Z <- MASS::mvrnorm(number_of_observations, mu = mu, Sigma = sigma_matrix)
#' S <- cov(Z) # Assume we have to estimate the mean
#'
#' g <- gips(S, number_of_observations)
#' if (require("graphics")) {
#'   plot(g, type = "MLE")
#' }
#'
#' g_map <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "hill_climbing")
#' if (require("graphics")) {
#'   plot(g_map, type = "both", logarithmic_x = TRUE)
#' }
#'
#' if (require("graphics")) {
#'   plot(g_map, type = "MLE")
#' }
#' # Now, the output is (most likely) different because the permutation
#'   # `g_map[[1]]` is (most likely) not an identity permutation.
#' 
#' g_map_MH <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "MH")
#' if (require("graphics")) {
#'   plot(g_map_MH, type = "n0")
#' }
plot.gips <- function(x, type = NA,
                      logarithmic_y = TRUE, logarithmic_x = FALSE,
                      color = NULL,
                      title_text = "Convergence plot",
                      xlabel = NULL, ylabel = NULL,
                      show_legend = TRUE,
                      ylim = NULL, xlim = NULL, ...) {
  if (!requireNamespace("graphics", quietly = TRUE)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "Package 'graphics' must be installed to use this function.",
      "x" = "Package 'graphics' seems to be unavailable."
    ))
  }

  validate_gips(x)

  if (length(type) != 1) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`type` must be an character vector of length 1.",
      "x" = paste0("You provided `type` with length ", length(type), " which is wrong!")
    ))
  }
  if (is.na(type)) {
    type <- ifelse(is.null(attr(x, "optimization_info")),
      "heatmap",
      "both"
    )

    rlang::inform(c("You used the default value of the 'type' argument in `plot()` for gips object.",
      "i" = paste0(
        "The `type` was automatically set to `type = '",
        type, "'`."
      )
    ))
  }

  if (type == "MLE") {
    type <- "heatmap"
  }

  if (!(type %in% c("heatmap", "block_heatmap", "all", "best", "both", "n0"))) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`type` must be one of: c('heatmap', 'MLE', 'block_heatmap', 'all', 'best', 'both', 'n0').",
      "x" = paste0("You provided `type == ", type, "`."),
      "i" = "Did You misspell the 'type' argument?"
    ))
  }

  if (type != "block_heatmap" && type != "heatmap" &&
    is.null(attr(x, "optimization_info"))) {
    rlang::abort(
      c(
        "There was a problem identified with the provided arguments:",
        "i" = "For non-optimized `gips` objects only the `type = 'heatmap', 'MLE' or 'block_heatmap'` can be used.",
        "x" = paste0(
          "You did not optimized `x` and provided `type = '",
          type, "'`."
        ),
        "i" = paste0(
          "Did You want to call `x <- find_MAP(g)` and then `plot(x, type = '",
          type, "')`?"
        ),
        "i" = "Did You want to use `type = 'heatmap'`?"
      )
    )
  }

  # dispatch to the appropriate internal plotting function
  if (type == "heatmap" || type == "block_heatmap") {
    return(plot_gips_heatmap(x, type = type))
  }
  if (type %in% c("all", "best", "both")) {
    plot_gips_convergence(
      x,
      type = type,
      logarithmic_y = logarithmic_y, logarithmic_x = logarithmic_x,
      color = color,
      title_text = title_text,
      xlabel = xlabel, ylabel = ylabel,
      show_legend = show_legend,
      ylim = ylim, xlim = xlim,
      ...
    )
  }
  if (type == "n0") {
    plot_gips_n0(
      x,
      logarithmic_y = logarithmic_y, logarithmic_x = logarithmic_x,
      color = color,
      title_text = title_text,
      xlabel = xlabel, ylabel = ylabel,
      show_legend = show_legend,
      ylim = ylim, xlim = xlim,
      ...
    )
  }

  invisible(NULL)
}


#' Plot a heatmap of the projected covariance matrix
#'
#' Internal function called by [plot.gips()] for
#' `type = "heatmap"` or `type = "block_heatmap"`.
#'
#' @param x A `gips` object.
#' @param type One of `"heatmap"` or `"block_heatmap"`.
#'
#' @returns A `ggplot` object.
#'
#' @noRd
plot_gips_heatmap <- function(x, type) {
  rlang::check_installed(c("dplyr", "tidyr", "tibble", "ggplot2"),
    reason = "to use `plot.gips()` with `type %in% c('heatmap', 'MLE', 'block_heatmap')`"
  )

  if (type == "block_heatmap") {
    my_projected_matrix <- get_diagonalized_matrix_for_heatmap(x)
  } else {
    my_projected_matrix <- project_matrix(attr(x, "S"), x[[1]])
  }

  p <- ncol(my_projected_matrix)

  if (is.null(colnames(my_projected_matrix))) {
    colnames(my_projected_matrix) <- as.character(seq_len(p))
  }
  if (is.null(rownames(my_projected_matrix))) {
    rownames(my_projected_matrix) <- as.character(seq_len(p))
  }

  my_rownames <- rownames(my_projected_matrix)
  my_colnames <- colnames(my_projected_matrix)
  rownames(my_projected_matrix) <- as.character(seq_len(p))
  colnames(my_projected_matrix) <- as.character(seq_len(p))

  # With this line, the R CMD check's "no visible binding for global variable" warning will not occur:
  col_id <- covariance <- row_id <- NULL

  # Life would be easier with pipes (%>%)
  my_transformed_matrix <- tibble::rownames_to_column(
    as.data.frame(my_projected_matrix),
    "row_id"
  )
  my_transformed_matrix <- tidyr::pivot_longer(my_transformed_matrix,
    -c(row_id),
    names_to = "col_id",
    values_to = "covariance"
  )
  my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
    col_id = as.numeric(col_id),
    row_id = as.numeric(row_id)
  )

  ggplot2::ggplot(
    my_transformed_matrix,
    ggplot2::aes(x = col_id, y = row_id, fill = covariance)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(na.value = "white") +
    ggplot2::scale_x_continuous(breaks = seq_len(p), labels = my_rownames) +
    ggplot2::scale_y_reverse(breaks = seq_len(p), labels = my_colnames) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Estimated covariance matrix\nprojected on permutation ", x[[1]]),
      x = "", y = ""
    )
}


#' Plot the optimizer posterior trajectory
#'
#' Internal function called by [plot.gips()] for
#' `type %in% c("all", "best", "both")`.
#'
#' @param x A `gips` object (must be optimized).
#' @param type One of `"all"`, `"best"`, or `"both"`.
#' @param logarithmic_y,logarithmic_x Booleans. Sets axes to log scale.
#' @param color Vector of line colors.
#' @param title_text Plot title.
#' @param xlabel,ylabel Axis labels.
#' @param show_legend Boolean. Whether to draw the legend.
#' @param ylim,xlim Axis limits; `NULL` uses data range.
#' @param ... Additional arguments passed to base graphics functions.
#'
#' @returns `NULL` invisibly.
#'
#' @noRd
plot_gips_convergence <- function(x, type,
                                  logarithmic_y, logarithmic_x,
                                  color,
                                  title_text,
                                  xlabel, ylabel,
                                  show_legend,
                                  ylim, xlim, ...) {
  if (is.null(ylabel)) {
    ylabel <- ifelse(logarithmic_y,
      "log posteriori",
      "posteriori"
    )
  }
  if (is.null(xlabel)) {
    xlabel <- ifelse(logarithmic_x,
      "log10 of number of function calls",
      "number of function calls"
    )
  }
  if (is.null(color)) {
    if (type == "both") {
      color <- c("red", "blue")
    } else {
      color <- "red"
    }
  }
  if (logarithmic_y) {
    y_values_from <- attr(x, "optimization_info")[["log_posteriori_values"]] # values of log_posteriori are logarithmic by default
  } else {
    y_values_from <- exp(attr(x, "optimization_info")[["log_posteriori_values"]])
  }

  y_values_max <- cummax(y_values_from)
  y_values_all <- y_values_from

  num_of_steps <- length(y_values_max)

  if (is.null(xlim)) {
    xlim <- c(1, num_of_steps)
  }

  if (is.null(ylim)) {
    ylim_plot <- c(min(y_values_from), y_values_max[num_of_steps])
    if (type == "best") {
      ylim_plot[1] <- y_values_from[1] # for the "best" type this is the smallest point of the graph
    }
  } else {
    ylim_plot <- ylim
  }

  # make the plot stairs-like
  x_points <- c(1, rep(2:num_of_steps, each = 2))

  if (logarithmic_x) {
    x_points <- log10(x_points)
    xlim <- log10(xlim)
  }

  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)

  if (type != "best") {
    # make the plot stairs-like
    y_points <- c(
      rep(y_values_all[1:(length(y_values_all) - 1)], each = 2),
      y_values_all[length(y_values_all)]
    )

    graphics::lines.default(x_points, y_points,
      type = "l", lwd = 3,
      col = color[1], # the first color
      ...
    )
  }
  if (type != "all") {
    # make the plot stairs-like
    y_points <- c(
      rep(y_values_max[1:(length(y_values_max) - 1)], each = 2),
      y_values_max[length(y_values_max)]
    )

    graphics::lines.default(x_points, y_points,
      lwd = 3, lty = 1,
      col = color[length(color)], # the last color
      ...
    )
  }

  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)

  if (show_legend) {
    if (type == "both") {
      legend_text <- c(
        "All calculated a posteriori",
        "Maximum a posteriori calculated"
      )
      lty <- c(1, 1)
      lwd <- c(3, 3)
    } else if (type == "all") {
      legend_text <- c("All calculated function values")
      lty <- 1
      lwd <- 3
    } else if (type == "best") {
      legend_text <- c("Maximum function values calculated")
      lty <- 1
      lwd <- 3
    }

    graphics::legend("bottomright",
      inset = .002,
      legend = legend_text,
      col = color,
      lty = lty, lwd = lwd,
      cex = 0.7, box.lty = 0
    )
  }

  invisible(NULL)
}


#' Plot the n0 trajectory from MH optimization
#'
#' Internal function called by [plot.gips()] for `type = "n0"`.
#'
#' @param x A `gips` object (must be optimized with MH).
#' @param logarithmic_y,logarithmic_x Booleans. Sets axes to log scale.
#' @param color Line color.
#' @param title_text Plot title.
#' @param xlabel,ylabel Axis labels.
#' @param show_legend Boolean. Whether to draw the legend.
#' @param ylim,xlim Axis limits; `NULL` uses data range.
#' @param ... Additional arguments passed to base graphics functions.
#'
#' @returns `NULL` invisibly.
#'
#' @noRd
plot_gips_n0 <- function(x,
                         logarithmic_y, logarithmic_x,
                         color,
                         title_text,
                         xlabel, ylabel,
                         show_legend,
                         ylim, xlim, ...) {
  if (is.null(ylabel)) {
    ylabel <- ifelse(logarithmic_y,
      "log n0",
      "n0"
    )
  }
  if (is.null(xlabel)) {
    xlabel <- ifelse(logarithmic_x,
      "log10 of number of function calls",
      "number of function calls"
    )
  }
  if (is.null(color)) {
    color <- "red"
  }

  if (logarithmic_y) {
    y_values <- log(attr(x, "optimization_info")[["all_n0"]])
  } else {
    y_values <- attr(x, "optimization_info")[["all_n0"]]
  }

  num_of_steps <- length(y_values)

  if (is.null(xlim)) {
    xlim <- c(1, num_of_steps)
  }

  if (is.null(ylim)) {
    ylim_plot <- c(0, max(y_values))
  } else {
    ylim_plot <- ylim
  }

  # make the plot stairs-like
  x_points <- c(1, rep(2:num_of_steps, each = 2))

  if (logarithmic_x) {
    x_points <- log10(x_points)
    xlim <- log10(xlim)
  }

  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)

  # make the plot stairs-like
  y_points <- c(
    rep(y_values[1:(length(y_values) - 1)], each = 2),
    y_values[length(y_values)]
  )

  graphics::lines.default(x_points, y_points,
    type = "l", lwd = 3,
    col = color[1], # the first color
    ...
  )

  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)

  if (show_legend) {
    legend_text <- "all perms n0"
    lty <- c(1, 1)
    lwd <- c(3, 3)

    graphics::legend("topright",
      inset = .002,
      legend = legend_text,
      col = color,
      lty = lty, lwd = lwd,
      cex = 0.7, box.lty = 0
    )
  }

  invisible(NULL)
}


#' Replace all non-block entries with NA
#'
#' Diagonalize matrix using found permutation and
#' replace all entries outside blocks (equal to 0) with NA.
#' This is done, because later these fields are plotted with background color.
#' It is more clear then.
#'
#' @param g `gips` object.
#' @noRd
get_diagonalized_matrix_for_heatmap <- function(g) {
  perm <- g[[1]]
  projected_matrix <- project_matrix(attr(g, "S"), perm)
  diagonalising_matrix <- prepare_orthogonal_matrix(perm)
  full_block_matrix <- t(diagonalising_matrix) %*% projected_matrix %*% diagonalising_matrix
  block_ends <- get_block_ends(get_structure_constants(perm))
  block_starts <- c(1, block_ends[-length(block_ends)] + 1)
  block_matrix <- matrix(
    nrow = nrow(full_block_matrix),
    ncol = ncol(full_block_matrix)
  )
  for (i in 1:length(block_starts)) {
    slice <- block_starts[i]:block_ends[i]
    block_matrix[slice, slice] <- full_block_matrix[slice, slice, drop = FALSE]
  }
  block_matrix
}
