#' Printing gips object
#' 
#' Printing function for gips class.
#' @param x object of class gips.
#' @param log_value logical. Weather to print the exp of a value of a \code{\link{goal_function}} or leave it in logarithmic form.
#' @param ... additional arguments passed to \code{\link{cat}}.
#' 
#' @return Invisible NULL.
#' @export
print.gips <- function(x, log_value = FALSE, ...){
  # TODO(Change it)
  
  value_part <- ifelse(log_value,
                       paste0(" with function log value ",
                              x[["found_point_function_logvalue"]]),
                       paste0(" with function value ",
                              exp(x[["found_point_function_logvalue"]])))
  cat(paste0("Optimization algorithm after ",
             length(x[["goal_function_logvalues"]]),
             " iterations found permutation ",
             x[["found_point"]],
             value_part),
      ...)
}




#' Plot optimization convergence
#' 
#' Plot method for gips objects.
#' @param x Object of class gips.
#' @param type Character. A type of a plot. Either "all", "best" or "both". For "all" plots values of goal function for all visited state. For "best" the best value of goal function up to the point are plotted. For "both" both lines are plotted.
#' @param logarithmic boolean.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ylim Limits of y axis. When \code{NULL}, the minimum and maximum of the \code{\link{goal_function}} is taken.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.gips <- function(x, type="both",
                      logarithmic=TRUE,
                      title_text="Convergence plot",
                      xlabel="number of function calls",
                      ylabel=NULL, show_legend=TRUE,
                      ylim=NULL, ...){
  stopifnot(type %in% c("all", "best", "both"))
  if(is.null(ylabel)){
    ylabel <- ifelse(logarithmic,
                     "log of a function",
                     "value of a function")
  }
  if(logarithmic){  # values of goal function are logarithmic by default
    y_values_from <- x[["goal_function_logvalues"]]
  }else{
    y_values_from <- exp(x[["goal_function_logvalues"]])
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
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)
  if(type != "best"){
    graphics::lines.default(1:num_of_steps, y_values_all, type = "l",
                            col = "blue", ...)
  }
  if(type != "all"){
    graphics::lines.default(1:num_of_steps, y_values_max, lwd=3,
                            lty = 2, col = "red", ...)
  }
  
  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)
  
  if(show_legend)
    graphics::legend("bottomright", inset=.002,
                     legend = c("All calculated function values",
                                "Maximum function values calculated"),
                     col = c("blue", "red"), lty = 1:2, cex = 0.7,
                     box.lty=0, lwd = c(1, 3))
  
  invisible(NULL)
}


# TODO(summary)


#' Plot optimization convergence for `best_growth` algorithm
#' 
#' Plot method for gips objects.
#' @param x Object of class gips.
#' @param logarithmic boolean.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ylim Limits of y axis. When \code{NULL}, the minimum and maximum of the \code{\link{goal_function}} is taken.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.optimized_best_growth <- function(x, logarithmic=TRUE,
                      title_text="Convergence plot",
                      xlabel="number of function calls",
                      ylabel=NULL, show_legend=TRUE,
                      ylim=NULL, ...){
  if(is.null(ylabel)){
    ylabel <- ifelse(logarithmic,
                     "log of a function",
                     "value of a function")
  }
  if(logarithmic){  # values of goal function are logarithmic by default
    y_values_from <- x[["goal_function_logvalues"]]
  }else{
    y_values_from <- exp(x[["goal_function_logvalues"]])
  }
  
  # for best_growth algorithm, x[["goal_function_logvalues"]] always increases
  y_values_all <- y_values_from
  
  num_of_neighbours <- choose(dim(x$U_used)[1], 2)
  num_of_steps <- (length(y_values_all) - 1) * num_of_neighbours + 1 # notice that when algorithm did not converged, x[["iterations_performed"]] + 1 == length(y_values_max); + 1 is for the starting point
  
  xlim <- c(1, num_of_steps)
  if(is.null(ylim))
    ylim_plot <- c(min(y_values_from), y_values_all[length(y_values_all)])
  else
    ylim_plot <- ylim
  
  x_axis_points <- (0:(length(y_values_all)-1)) * num_of_neighbours + 1
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim_plot)
  graphics::lines.default(x_axis_points, y_values_all, type = "l",
                          col = "blue", ...)
  
  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)
  
  if(show_legend)
    graphics::legend("bottomright", inset=.002,
                     legend = c("Function values of visited nodes"),
                     col = "blue", lty = 1, cex = 0.7,
                     box.lty=0, lwd = 1)
  
  invisible(NULL)
}




