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
#'     
#' Arguments `logarithmic_y`, `logarithmic_x`, `color`, `title_text`,
#'     `xlabel`, `ylabel`, `show_legend`, `ylim`, and `xlim` are only used for
#'     `type %in% c("all", "best", "both", "n0")` and ignored for heatmap types.
#' @param logarithmic_y,logarithmic_x A boolean.
#'     Sets the axis of the plot in logarithmic scale.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param color Vector of colors to be used to plot lines.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param title_text Text to be in the title of the plot.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param xlabel Text to be on the bottom of the plot.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param ylabel Text to be on the left of the plot.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param show_legend A boolean. Whether or not to show a legend.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param ylim Limits of the y axis. When `NULL`, uses the data range.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param xlim Limits of the x axis. When `NULL`, uses the data range.
#'     Only used for `type %in% c("all", "best", "both", "n0")`.
#' @param ... Ignored.
#'
#' @returns An object of class `ggplot`.
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
#' plot(g, type = "MLE")
#'
#' g_map <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "hill_climbing")
#' plot(g_map, type = "both", logarithmic_x = TRUE)
#'
#' plot(g_map, type = "MLE")
#' # Now, the output is (most likely) different because the permutation
#'   # `g_map[[1]]` is (most likely) not an identity permutation.
#'
#' g_map_MH <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "MH")
#' plot(g_map_MH, type = "n0", logarithmic_y = FALSE)
plot.gips <- function(x, type = NA,
                      logarithmic_y = TRUE, logarithmic_x = FALSE,
                      color = NULL,
                      title_text = "Convergence plot",
                      xlabel = NULL, ylabel = NULL,
                      show_legend = TRUE,
                      ylim = NULL, xlim = NULL, ...) {
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

  # Check for required plotting packages
  rlang::check_installed(c("ggplot2", "dplyr", "tidyr", "tibble"),
    reason = "to use `plot.gips()`"
  )

  # dispatch to the appropriate internal plotting function
  if (type %in% c("heatmap", "block_heatmap")) {
    return(plot_gips_heatmap(x, type = type))
  }
  if (type %in% c("all", "best", "both")) {
    return(plot_gips_convergence(
      x,
      type = type,
      logarithmic_y = logarithmic_y, logarithmic_x = logarithmic_x,
      color = color,
      title_text = title_text,
      xlabel = xlabel, ylabel = ylabel,
      show_legend = show_legend,
      ylim = ylim, xlim = xlim
    ))
  }
  if (type == "n0") {
    # Check that optimization was done with MH
    opt_algo <- attr(x, "optimization_info")$optimization_algorithm_used
    is_mh <- if (length(opt_algo) == 1) {
      opt_algo == "Metropolis_Hastings"
    } else {
      opt_algo[-1] == "Metropolis_Hastings"
    }
    
    if (!is_mh) {
      rlang::abort(c(
        "There was a problem identified with the provided arguments:",
        "i" = "`type = 'n0'` can only be used with `gips` objects optimized using the Metropolis-Hastings algorithm.",
        "x" = paste0(
          "Your `gips` object was optimized using: ",
          paste(opt_algo, collapse = " -> "), "."
        ),
        "i" = "Did you mean to use `type = 'all'`, `type = 'best'`, or `type = 'both'` instead?"
      ))
    }
    
    return(plot_gips_n0(
      x,
      logarithmic_y = logarithmic_y, logarithmic_x = logarithmic_x,
      color = color,
      title_text = title_text,
      xlabel = xlabel, ylabel = ylabel,
      show_legend = show_legend,
      ylim = ylim, xlim = xlim
    ))
  }
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
#'
#' @returns A `ggplot` object.
#'
#' @noRd
plot_gips_convergence <- function(x, type,
                                  logarithmic_y, logarithmic_x,
                                  color,
                                  title_text,
                                  xlabel, ylabel,
                                  show_legend,
                                  ylim, xlim) {
  if (is.null(ylabel)) {
    ylabel <- ifelse(logarithmic_y, "log posteriori", "posteriori")
  }
  if (is.null(xlabel)) {
    xlabel <- ifelse(logarithmic_x,
      "log10 of number of function calls",
      "number of function calls"
    )
  }

  log_post <- attr(x, "optimization_info")[["log_posteriori_values"]] # stored in log scale
  y_values_from <- if (logarithmic_y) log_post else exp(log_post)
  y_values_max  <- cummax(y_values_from)
  num_of_steps  <- length(y_values_from)

  all_label  <- "All calculated a posteriori"
  best_label <- "Maximum a posteriori calculated"

  # R CMD check: no visible binding for global variable
  step <- value <- series <- NULL

  df_all  <- data.frame(step = seq_len(num_of_steps), value = y_values_from, series = all_label)
  df_best <- data.frame(step = seq_len(num_of_steps), value = y_values_max,  series = best_label)

  df <- switch(type,
    "all"  = df_all,
    "best" = df_best,
    "both" = rbind(df_all, df_best)
  )

  if (is.null(color)) {
    color <- if (type == "both") c("red", "blue") else "red"
  }
  color_map <- switch(type,
    "all"  = stats::setNames(color[1],              all_label),
    "best" = stats::setNames(color[1],              best_label),
    "both" = stats::setNames(c(color[1], color[2]), c(all_label, best_label))
  )

  if (is.null(xlim)) xlim <- c(1, num_of_steps)
  if (is.null(ylim)) {
    y_min <- if (type == "best") y_values_from[1] else min(y_values_from)
    ylim  <- c(y_min, y_values_max[num_of_steps])
  }

  g_plot <- ggplot2::ggplot(df, ggplot2::aes(x = step, y = value, color = series)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = title_text, x = xlabel, y = ylabel, color = NULL) +
    ggplot2::theme_bw()

  if (logarithmic_x) {
    g_plot <- g_plot + ggplot2::scale_x_log10()
  }

  if (!show_legend) {
    g_plot <- g_plot + ggplot2::theme(legend.position = "none")
  }

  g_plot
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
#'
#' @returns A `ggplot` object.
#'
#' @noRd
plot_gips_n0 <- function(x,
                         logarithmic_y, logarithmic_x,
                         color,
                         title_text,
                         xlabel, ylabel,
                         show_legend,
                         ylim, xlim) {
  if (is.null(ylabel)) {
    ylabel <- ifelse(logarithmic_y, "log n0", "n0")
  }
  if (is.null(xlabel)) {
    xlabel <- ifelse(logarithmic_x,
      "log10 of number of function calls",
      "number of function calls"
    )
  }
  if (is.null(color)) color <- "red"

  raw_n0 <- attr(x, "optimization_info")[["all_n0"]]
  y_values <- if (logarithmic_y) log(raw_n0) else raw_n0
  num_of_steps <- length(y_values)

  if (is.null(xlim)) xlim <- c(1, num_of_steps)
  if (is.null(ylim)) ylim <- c(1, max(y_values))

  n0_label <- "all perms n0"

  # R CMD check: no visible binding for global variable
  step <- value <- label <- NULL

  df <- data.frame(step = seq_len(num_of_steps), value = y_values, label = n0_label)

  g_plot <- ggplot2::ggplot(df, ggplot2::aes(x = step, y = value, color = label)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::scale_color_manual(values = stats::setNames(color[1], n0_label)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = title_text, x = xlabel, y = ylabel, color = NULL) +
    ggplot2::theme_bw()

  if (logarithmic_x) {
    g_plot <- g_plot + ggplot2::scale_x_log10()
  }

  if (!show_legend) {
    g_plot <- g_plot + ggplot2::theme(legend.position = "none")
  }

  g_plot
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
