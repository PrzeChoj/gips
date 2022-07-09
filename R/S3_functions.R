#' Printing gips object
#' 
#' Printing function for gips class.
#' @param x object of class gips.
#' @param log_value logical. Weather to print the exp of a value of a \code{\link{log_likelihood_of_perm}} or leave it in logarithmic form.
#' @param ... additional arguments passed to \code{\link{cat}}.
#' 
#' @return Invisible NULL.
#' @export
print.gips <- function(x, log_value = TRUE, ...){
  # TODO(it is not likelihood, but sth proportional to it. See #ISSUE11)
  
  value_part <- ifelse(log_value,
                       paste0(" with log likelihood ",
                              x[["found_perm_log_likelihood"]]),
                       paste0(" with likelihood ",
                              exp(x[["found_perm_log_likelihood"]])))
  cat(paste0("Optimization algorithm ",
             x[["optimization_algorithm_used"]], " after ",
             length(x[["log_likelihood_values"]]),
             " iterations found permutation ",
             x[["found_perm"]],
             value_part),
      ...)
}




#' Plot optimization convergence
#' 
#' Plot method for gips objects.
#' @param x Object of class gips.
#' @param type Character. A type of a plot. Either "all", "best" or "both". For "all" plots likelihood for all visited state. For "best" the biggest likelihood up to the point are plotted. For "both" both lines are plotted.
#' @param logarithmic_y boolean.
#' @param logarithmic_x boolean.
#' @param color Vector of olors to be used to plot lines.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ylim Limits of y axis. When \code{NULL}, the minimum and maximum of the \code{\link{log_likelihood_of_perm}} is taken.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.gips <- function(x, type="both",
                      logarithmic_y=TRUE,
                      logarithmic_x=FALSE,
                      color=NULL,
                      title_text="Convergence plot",
                      xlabel=NULL,
                      ylabel=NULL, show_legend=TRUE,
                      ylim=NULL, ...){
  # TODO(it is not likelihood, but sth proportional to it. See #ISSUE11)
  
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop(
      "Package \"graphics\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  stopifnot(type %in% c("all", "best", "both"))
  if(is.null(ylabel)){
    ylabel <- ifelse(logarithmic_y,
                     "log likelihood",
                     "likelihood")
  }
  if(is.null(xlabel)){
    xlabel <- ifelse(logarithmic_x,
                     "log10 of number of function calls",
                     "number of function calls")
  }
  if(is.null(color)){
    if(type == "both"){
      color <- c("blue", "red")
    }else{
      color <- "blue"
    }
  }
  if(logarithmic_y){
    y_values_from <- x[["log_likelihood_values"]] # values of likelihood are logarithmic by default
  }else{
    y_values_from <- exp(x[["log_likelihood_values"]])
  }

  y_values_max <- cummax(y_values_from)
  y_values_all <- y_values_from
  
  num_of_steps <- length(y_values_max)
  
  xlim <- c(1, num_of_steps)
  if(is.null(ylim)){
    ylim_plot <- c(min(y_values_from), y_values_max[num_of_steps])
    if(type == "best"){
      ylim_plot[1] <- y_values_from[1] # for the "best" type this is the smallest point of the graph
    }
  }
  else
    ylim_plot <- ylim
  
  # make the plot stairs-like
  x_points <- c(1, rep(2:num_of_steps, each = 2))
  
  if(logarithmic_x){
    x_points <- log10(x_points)
    xlim <- log10(xlim)
  }
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)
  
  if(type != "best"){
    # make the plot stairs-like
    y_points <- c(rep(y_values_all[1:(length(y_values_all)-1)], each = 2),
                  y_values_all[length(y_values_all)])
    
    graphics::lines.default(x_points, y_points, type = "l",
                            lwd=3, col = color[1], # the first color
                            ...)
  }
  if(type != "all"){
    # make the plot stairs-like
    y_points <- c(rep(y_values_max[1:(length(y_values_max)-1)], each = 2),
                  y_values_max[length(y_values_max)])
    
    graphics::lines.default(x_points, y_points, lwd=3,
                            lty = 1,
                            col = color[length(color)], # the last color
                            ...)
  }
  
  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)
  
  if(show_legend){
    if(type == "both"){
      legend_text <- c("All calculated likelihoods",
                       "Maximum likelihoods calculated")
      lty <- c(1, 1)
      lwd <- c(3, 3)
    }else if(type == "all"){
      legend_text <- c("All calculated function values")
      lty <- 1
      lwd <- 3
    }else if(type == "best"){
      legend_text <- c("Maximum function values calculated")
      lty <- 1
      lwd <- 3
    }
    
    graphics::legend("bottomright", inset=.002,
                     legend = legend_text,
                     col = color, lty = lty, cex = 0.7,
                     box.lty=0, lwd = lwd)
  }
  
  invisible(NULL)
}


# TODO(summary, the n_0 of best perms, the distribution of likelihood values)
summary.gips <- function(object, ...){
  invisible(NULL)
}




