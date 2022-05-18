#' Printing MH
#' 
#' Printing function for MH.
#' @param x object of class MH. Output of MH function.
#' @param ... additional arguments passed to \code{\link{cat}}.
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
#' @param x object of class MH. Output of MH function.
#' @param logarithmic boolean.
#' @param type To plot the maximum found value or all found values.
#' @param title_text Text to be in a title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param ... additional arguments passed to \code{\link{print}}.
#' 
#' @return Invisible NULL.
#' @export
plot.MH <- function(x, logarithmic=TRUE, type="max",
                    title_text="Convergence plot",
                    xlabel="number of function calls",
                    ylabel=NULL, ...){
  stopifnot("MH" %in% class(x))
  if(is.null(ylabel)){
    ylabel <- ifelse(logarithmic,
                     "biggest log of a function",
                     "biggest value of a function")
  }
  if(logarithmic){  # values of goal function are logarithmic by default
    y_values_from <- x[["goal_function_logvalues"]]
  }else{
    y_values_from <- exp(x[["goal_function_logvalues"]])
  }
  if(type == "max"){
    y_values <- cummax(y_values_from)
  }else{
    y_values <- y_values_from
  }
  
  xlim <- c(1, length(x[["goal_function_logvalues"]]))
  ylim <- range(y_values)
  
  graphics::plot.new()
  graphics::plot.window(xlim, ylim)
  graphics::lines.default(1:length(y_values), y_values, ...)
  
  graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
  graphics::axis(1, ...)
  graphics::axis(2, ...)
  graphics::box(...)
  
  invisible(NULL)
}



