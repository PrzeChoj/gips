#' Printing MH
#' 
#' Printing function for MH.
#' @param x object of class MH. Output of MH function.
#' @param ... additional arguments passed to \code{\link{cat}}.
#' 
#' @return Invisible NULL.
#' @export
print.MH <- function(x, ...){
  # TODO(Change it)
  cat(paste0("Metropolis-Hastings algorithm after ",
             length(x[["goal_function_logvalues"]]),
             " iterations found permutation ",
             x[["found_point"]],
             " with function value ",
             exp(x[["found_point_function_logvalue"]])),
      ...)
}




#' Plot Metropolis-Hasting convergence
#' 
#' Plot method for MH objects returned by \code{\link{MH}}.
#' @param x Object of class MH. Output of MH function.
#' @param logarithmic boolean.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend boolean.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.MH <- function(x, logarithmic=TRUE,
                    title_text="Convergence plot",
                    xlabel="number of function calls",
                    ylabel=NULL, show_legend=TRUE,
                    ...){
  stopifnot("MH" %in% class(x))
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
  ylim <- c(y_values_max[1], y_values_max[num_of_steps])
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim)
  graphics::lines.default(1:num_of_steps, y_values_all, type = "l",
                          col = "blue", ...)
  graphics::lines.default(1:num_of_steps, y_values_max, lwd=3,
                          lty = 2, col = "red", ...)
  
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



